#include "runAction.hh"
#include "construction.hh"
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <sstream>
#include <ctime>
#include <cmath>
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4Threading.hh"

namespace
{
	inline void ensure_dir(const std::string& path)
	{
		std::error_code ec;
		if (!std::filesystem::exists(path)) 
		{
			std::filesystem::create_directories(path, ec);
		}
		if (ec) 
		{
			G4cerr << "[RunAction] Failed to create dir: " << path << " : " << ec.message() << G4endl;
		}
  	}

  	// write 1D vector as CSV 
	template<typename T>
	void write_csv_1d(const std::string& path, const std::vector<T>& v, double scale = 1.0) 
	{
		if (v.empty()) return;
		std::ofstream ofs(path);
		ofs.setf(std::ios::scientific, std::ios::floatfield);
		ofs << std::setprecision(3);
		for (std::size_t i = 0; i < v.size(); ++i)
		{
			if (i) ofs << ',';
			ofs << (static_cast<double>(v[i]) * scale);
		}
		ofs << "\n";
	}
	// write 2D vector as CSV nx rows, nz columns
	// data is in row-major order: data[iz * nx + ix]
	// scale is a multiplicative factor applied to each element (e.g. for unit conversion)
	template<typename T>
	void write_csv_2d(const std::string& path, const std::vector<T>& data, int nx, int nz, double scale = 1.0)
	{
    	if (data.size() != static_cast<size_t>(nx) * static_cast<size_t>(nz) || nx <= 0 || nz <= 0) return;

    	std::ofstream ofs(path);
    	ofs.setf(std::ios::scientific, std::ios::floatfield);
    	ofs << std::setprecision(3);

    // Outer: x (rows). Inner: z (columns, left->right).
    for (int ix = 0; ix < nx; ++ix) {
        for (int iz = 0; iz < nz; ++iz) {
            if (iz) ofs << ',';
            ofs << (static_cast<double>(data[iz * nx + ix]) * scale);
        }
        ofs << '\n';
    }
}
	// YYYY-MM-DD (UTC)
	inline std::string utc_date() 
	{
		std::time_t t = std::time(nullptr);
		std::tm g{}; gmtime_r(&t, &g);
		char buf[16];
		std::strftime(buf, sizeof(buf), "%Y-%m-%d", &g);
		return buf;
    }
	// HH-MM-SS (UTC)
	inline std::string utc_time() 
	{
		std::time_t t = std::time(nullptr);
		std::tm g{}; gmtime_r(&t, &g);
		char buf[16];
		std::strftime(buf, sizeof(buf), "%H-%M-%S", &g);
		return buf;
    }	
	// Safe file copy
  	inline void try_copy_file(const std::string& src, const std::string& dst) 
	{
		if (src.empty()) return;
		std::error_code ec;
		if (!std::filesystem::exists(src)) 
		{
			G4cerr << "[RunAction] Source file not found (skip copy): " << src << G4endl;
			return;
		}
		std::filesystem::copy_file(src, dst,std::filesystem::copy_options::overwrite_existing, ec);
		if (ec) 
		{
			G4cerr << "[RunAction] Copy failed: " << src << " -> " << dst << " : " << ec.message() << G4endl;
		}
    }
	// Compute normalization factor for Maxwell-like spectrum with low-energy cutoff
	inline G4double ComputeSpectrumNormFactor(const SimConfig& cfg)
	{
		if (cfg.EminCutoff <= 0.)
			return 1.0;

		const G4double a  = cfg.EminCutoff / cfg.Tmaxwell; // dimensionless
		const G4double s  = std::sqrt(a);

		// F_below = erf(s) - 2/sqrt(pi)*s*exp(-a)
		const G4double F_below = std::erf(s) - (2.0 / std::sqrt(CLHEP::pi)) * s * std::exp(-a);

		G4double frac_above = 1.0 - F_below;
		if (frac_above <= 0.0) 
			return 1.0; // extreme cutoff, avoid blow-up

		return 1.0 / frac_above;
	}
}

MyRunAction::~MyRunAction()
{}

G4Run* MyRunAction::GenerateRun()
{
	fRun = new Run(cfg_);
	return fRun;

}

void MyRunAction::BeginOfRunAction(const G4Run* run)
{
	// Get angle ranges from detector construction (geometry is now built)
	if (cfg_.enableExitPlane) {
		auto* det = dynamic_cast<const MyDetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
		if (det) {
			fAngleMinThetaX = det->GetAngleMinThetaX();
			fAngleMaxThetaX = det->GetAngleMaxThetaX();
			fAngleMinThetaY = det->GetAngleMinThetaY();
			fAngleMaxThetaY = det->GetAngleMaxThetaY();
			if (IsMaster()) {
				G4cout << "[RunAction::Begin] Angle ranges: theta_x=[" << fAngleMinThetaX << "," << fAngleMaxThetaX 
				       << "], theta_y=[" << fAngleMinThetaY << "," << fAngleMaxThetaY << "]" << G4endl;
			}
		}
	}
	
	// Make folder: data/<YYYY-MM-DD>/<HH-MM-SS>+<runName>
	date_ = utc_date();
	timeStamp_  = utc_time();
	runName_ = cfg_.runName.empty() ? "run" : cfg_.runName;
	dateDir_ = cfg_.simOutDir + "/" + date_;
	ensure_dir(cfg_.simOutDir);
	ensure_dir(dateDir_);

	stamp_  = timeStamp_ + "__" + runName_;
	outDir_ = dateDir_ + "/" + stamp_;
	ensure_dir(outDir_);

	// Start timer
    myTimer_.Start(); 
	G4cout << "[RunAction] Output dir: " << outDir_ << G4endl;
	// Compute spectrum normalization factor
	spectrumNorm_ = ComputeSpectrumNormFactor(cfg_);
}
void MyRunAction::EndOfRunAction(const G4Run* run)
{
    if (!IsMaster()) return;
	 const G4int N_sim = run->GetNumberOfEvent();     // actually simulated
    const G4double N_phys = N_sim * spectrumNorm_;    // effective ones
	const Run* simRun = static_cast<const Run*>(run);
	// Stop timer
	myTimer_.Stop();
	const double dt = myTimer_.GetRealElapsed();
    
	try_copy_file(cfg_.cfgDir + "/" + cfg_.cfgFile , outDir_ + "/" + cfg_.cfgFile); // copy config file to output dir
	try_copy_file(macroPath_, outDir_ + "/macro_used.mac"); // copy macro to output dir
	
    // normalization for fluence maps (path length / area), and optionally per-primary
  	double scaleFluence = 1.0;
	double scaleFluenceWorld = 1.0;
  	if (cfg_.normalizeByArea)
    {
    	const double area_mm2 = (cfg_.pixelX / mm) * (cfg_.pixelZ / mm); // numeric mm^2
    	if (area_mm2 > 0.0) scaleFluence /= area_mm2;                    // length(mm)/area(mm^2) => 1/mm
		if (cfg_.runWorldMap)
		{
			const double areaWorld_mm2 = (cfg_.worldPixelX / mm) * (cfg_.worldPixelZ / mm); // numeric mm^2
			if (areaWorld_mm2 > 0.0) scaleFluenceWorld /= areaWorld_mm2;                    // length(mm)/area(mm^2) => 1/mm
		}
	}
  	if (cfg_.normalizePerPrimary)
	{
    	if (N_phys > 0) scaleFluence /= static_cast<double>(N_phys);
		if (cfg_.runWorldMap)
		{
			if (N_phys > 0) scaleFluenceWorld /= static_cast<double>(N_phys);
		}
  	}
	if (cfg_.runPiPlusMain)
	{
		 // spectra
    	write_csv_1d(outDir_ + "/pion_energy_in.csv",  simRun->pionEnergyIn, spectrumNorm_);
    	write_csv_1d(outDir_ + "/pion_energy_out.csv", simRun->pionEnergyOut,spectrumNorm_);

		// absorber-local pion fluence
		write_csv_2d(outDir_ + "/pion_fluence_abs.csv",
					simRun->pionFluenceAbs, cfg_.nAbsorberX, cfg_.nAbsorberZ, scaleFluence);

		// optional world-wide map (write only if the vector is non-empty) ADD OPTION TO CFG
		if (!simRun->pionFluenceWorld.empty())
		write_csv_2d(outDir_ + "/pion_fluence_world.csv", simRun->pionFluenceWorld, cfg_.nWorldX, cfg_.nWorldZ, scaleFluenceWorld);
	}
	// Converter stats (electrons/gammas)
 	if (cfg_.runConvStats)
   	{
		// gamma spectrum inside absorber (your “gamma_energy_in”)
		write_csv_1d(outDir_ + "/gamma_energy_in.csv", simRun->gammaEnergy, spectrumNorm_);

		// absorber-local fluence maps for e- and gamma
		write_csv_2d(outDir_ + "/e_fluence_abs.csv", simRun->eFluenceAbs, cfg_.nAbsorberX, cfg_.nAbsorberZ, scaleFluence);

		write_csv_2d(outDir_ + "/gamma_fluence_abs.csv", simRun->gammaFluenceAbs, cfg_.nAbsorberX, cfg_.nAbsorberZ, scaleFluence);

		write_csv_2d(outDir_ + "/gamma_fluence_over200_abs.csv", simRun->gammaFluenceOver200Abs, cfg_.nAbsorberX, cfg_.nAbsorberZ, scaleFluence);

		write_csv_2d(outDir_ + "/gamma_creation_abs.csv", simRun->gammaCreationAbs, cfg_.nAbsorberX, cfg_.nAbsorberZ, scaleFluence);
  	}
	if (cfg_.runDebug)
	{
		// genrator beam distribution in the X-Y plane of the absorber
		write_csv_2d(outDir_ + "/genrator_beam_xy.csv", simRun->genratorBeamXY, cfg_.nAbsorberX, cfg_.nAbsorberY);
		// genrator energy histogram
		write_csv_1d(outDir_ + "/genrator_energy.csv", simRun->genratorEnergy);
	}
	
	// Exit plane angle histograms (theta_x vs theta_y)
	if (cfg_.enableExitPlane)
	{
		write_csv_2d(outDir_ + "/pion_exit_plane_angles.csv", 
					 simRun->pionExitPlaneAngleHistograms, 
					 cfg_.nAngleBinsThetaX, cfg_.nAngleBinsThetaY);
		if (cfg_.angleIncludeBackground)
		{
			write_csv_2d(outDir_ + "/background_exit_plane_angles.csv", 
						 simRun->backgroundExitPlaneAngleHistograms, 
						 cfg_.nAngleBinsThetaX, cfg_.nAngleBinsThetaY);
		}
		
		// Write angle bin edges metadata
		std::ofstream angleMeta(outDir_ + "/angle_bin_edges.txt");
		angleMeta.setf(std::ios::scientific, std::ios::floatfield);
		angleMeta << std::setprecision(6);
		angleMeta << "# Exit plane angle histogram metadata\n";
		angleMeta << "# Angles are in radians, using atan2(px,pz) and atan2(py,pz)\n";
		angleMeta << "# Histogram dimensions: " << cfg_.nAngleBinsThetaX << " x " << cfg_.nAngleBinsThetaY << "\n\n";
		
		angleMeta << "[ThetaX Bins]\n";
		angleMeta << "min = " << fAngleMinThetaX << " rad\n";
		angleMeta << "max = " << fAngleMaxThetaX << " rad\n";
		angleMeta << "nbins = " << cfg_.nAngleBinsThetaX << "\n";
		angleMeta << "bin_width = " << ((fAngleMaxThetaX - fAngleMinThetaX) / cfg_.nAngleBinsThetaX) << " rad\n\n";
		
		angleMeta << "[ThetaY Bins]\n";
		angleMeta << "min = " << fAngleMinThetaY << " rad\n";
		angleMeta << "max = " << fAngleMaxThetaY << " rad\n";
		angleMeta << "nbins = " << cfg_.nAngleBinsThetaY << "\n";
		angleMeta << "bin_width = " << ((fAngleMaxThetaY - fAngleMinThetaY) / cfg_.nAngleBinsThetaY) << " rad\n\n";
		
		// Write bin centers for easy plotting
		angleMeta << "[ThetaX Bin Centers (rad)]\n";
		const double dThetaX = (fAngleMaxThetaX - fAngleMinThetaX) / cfg_.nAngleBinsThetaX;
		for (int i = 0; i < cfg_.nAngleBinsThetaX; ++i) 
		{
			angleMeta << (fAngleMinThetaX + (i + 0.5) * dThetaX);
			if (i < cfg_.nAngleBinsThetaX - 1) angleMeta << ", ";
		}
		angleMeta << "\n\n";
		
		angleMeta << "[ThetaY Bin Centers (rad)]\n";
		const double dThetaY = (fAngleMaxThetaY - fAngleMinThetaY) / cfg_.nAngleBinsThetaY;
		for (int i = 0; i < cfg_.nAngleBinsThetaY; ++i) 
		{
			angleMeta << (fAngleMinThetaY + (i + 0.5) * dThetaY);
			if (i < cfg_.nAngleBinsThetaY - 1) angleMeta << ", ";
		}
		angleMeta << "\n";
	}
	
	// Write results.txt summary

	std::ofstream txt(outDir_ + "/results.txt");
    txt.setf(std::ios::scientific, std::ios::floatfield);
    txt << std::setprecision(1);

    // Header
    txt << "Run folder: " << stamp_ << "\n";
    txt << "Date      : " << date_ << "\n\n";
	txt << "seed used : " << G4Random::getTheSeed() << "\n\n";

    // From cfg
    txt << "[Config]\n";
    txt << "  run name           : " << runName_ << "\n";
    txt << "  run Description    : " << cfg_.runDescription << "\n";
    txt.setf(std::ios::fixed, std::ios::floatfield);
    txt << std::setprecision(2);
    txt << "  world (half, cm)   : "
        << (cfg_.worldX/cm) << " x " << (cfg_.worldY/cm) << " x " << (cfg_.worldZ/cm) << "\n";
    txt << "  absorber(half, cm) : "
        << (cfg_.absorberX/cm) << " x " << (cfg_.absorberY/cm) << " x " << (cfg_.absorberZ/cm) << "\n";
    txt << "  magnetic field     : " << (cfg_.enableMagneticField ? "enabled" : "disabled") << "\n";
    if (cfg_.enableMagneticField) 
	{
        txt << "  field strength (T) : " << (cfg_.fieldZ/tesla) << "\n";
    }
    // TODO: field location (add to SimConfig later)
    // txt << "  field location    : " << cfg_.fieldLocation << "\n";
    txt << "  pixel (mm)         : " << (cfg_.pixelX/mm) << " x " << (cfg_.pixelZ/mm) << "\n";
    txt << "  primary particle   : " << cfg_.gunParticle << "\n";
	if (cfg_.beamModel == SimConfig::BeamModel::TopHatDisk)
	{
		txt << "  beam model         : TopHatDisk (radius = " << (cfg_.beamRadius/mm) << " mm)\n";
	}
	else if (cfg_.beamModel == SimConfig::BeamModel::Gaussian)
	{
		txt << "  beam model         : Gaussian (sigma = " << (cfg_.beamSigma/mm) << " mm)\n";
	}
	else
	{	
		txt << "  beam model         : Unknown\n";
	}
	if (cfg_.energyModel == SimConfig::EnergyModel::Mono)
	{
		txt << "  energy model       : Mono (E = " << (cfg_.Emono/MeV) << " MeV)\n";
	}
	else if (cfg_.energyModel == SimConfig::EnergyModel::Uniform)
	{
		txt << "  energy model       : Uniform (E = " << (cfg_.EunifMin/MeV) << " .. " << (cfg_.Emax/MeV) << " MeV)\n";
	}
	else if (cfg_.energyModel == SimConfig::EnergyModel::MaxwellLike)
	{
		txt << "  energy model       : Maxwell-like (T = " << (cfg_.Tmaxwell/MeV) << " MeV, Emin cutoff = " << (cfg_.EminCutoff/MeV) << " MeV)\n";
	}
	else
	{
		txt << "  energy model       : Unknown\n";
	}
	if (cfg_.enableExitPlane)
	{
		txt << "  exit plane Z (cm)  : " << (cfg_.exitPlaneZ/cm) << "\n";
		txt << "  exit plane size    : " << (cfg_.exitPlaneHalfX/cm) << " x " << (cfg_.exitPlaneHalfY/cm) << " cm\n";
		txt << "  angle bins         : " << cfg_.nAngleBinsThetaX << " x " << cfg_.nAngleBinsThetaY << "\n";
		txt.setf(std::ios::fixed, std::ios::floatfield);
		txt << std::setprecision(4);
		txt << "  theta_x range (rad): " << fAngleMinThetaX << " to " << fAngleMaxThetaX 
		    << " (" << (fAngleMinThetaX * 180.0 / M_PI) << "° to " << (fAngleMaxThetaX * 180.0 / M_PI) << "°)\n";
		txt << "  theta_y range (rad): " << fAngleMinThetaY << " to " << fAngleMaxThetaY 
		    << " (" << (fAngleMinThetaY * 180.0 / M_PI) << "° to " << (fAngleMaxThetaY * 180.0 / M_PI) << "°)\n";
		txt << "  angle auto-compute : " << (cfg_.angleAutoCompute ? "yes" : "no") << "\n";
		txt.setf(std::ios::scientific, std::ios::floatfield);
		txt << std::setprecision(1);
	}

    txt << "\n";

    // From macro (runtime parameters)
    txt << "[Macro]\n";
    if (N_sim >= 10000)
	{
        txt.setf(std::ios::scientific, std::ios::floatfield);
        txt << std::setprecision(1);
        txt << "  number of primaries      : " << static_cast<double>(N_sim) << "\n";
        txt << "  effective primaries (Emin=" << (cfg_.EminCutoff/MeV) << " MeV): " << N_phys << "\n";
    } 
	else
	{
        txt.setf(std::ios::fixed, std::ios::floatfield);
        txt << std::setprecision(0);
        txt << "  number of primaries      : " << N_sim << "\n";
        txt << std::setprecision(1);
        txt << "  effective primaries (Emin=" << (cfg_.EminCutoff/MeV) << " MeV): " << N_phys << "\n";
    }
    txt.setf(std::ios::fixed, std::ios::floatfield);
    txt << std::setprecision(0);
    txt << "  number of threads        : " << G4RunManager::GetRunManager()->GetNumberOfThreads() << "\n";
    if (!macroPath_.empty())
	{
    	txt << "  macro file         : " << macroPath_ << "\n";
	}
    txt << "\n";

    // From simulation counters
    txt << "[Results]\n";
    txt.setf(std::ios::fixed, std::ios::floatfield);
    txt << std::setprecision(0);
    txt << "  pi+ created (in absorber): " << simRun->npiPosIn << "\n";
    txt << "  pi+ exited absorber      : " << simRun->npiPosOut << "\n";
    txt << "\n";

    // Timing
    const auto nEvt = static_cast<double>(N_sim);
    const double perEvt = (nEvt > 0.0) ? (dt / nEvt) : 0.0;
    txt << "[Timing]\n";
    txt.setf(std::ios::fixed, std::ios::floatfield);
    txt << std::setprecision(2);
    txt << "  runtime (s)          : " << dt << "\n";
    txt << "  runtime per primary (s): ";
    if (perEvt < 0.0001)
	{
        txt.setf(std::ios::scientific, std::ios::floatfield);
        txt << std::setprecision(2) << perEvt << "\n";
        txt.setf(std::ios::fixed, std::ios::floatfield);
    } 
	else
	{
        txt << perEvt << "\n";
    }
    txt << "  system time (s)      : " << myTimer_.GetSystemElapsed() << "\n";
    txt << "  user time (s)        : " << myTimer_.GetUserElapsed() << "\n";
    txt << "  cpu efficiency (%)   : " << (100.0 * (myTimer_.GetUserElapsed()+myTimer_.GetSystemElapsed() )/ dt) << "\n";
    txt << "\n";



	G4cout
    << "\n----- End of Run (master) -----\n"
	<< "  Output dir       : " << outDir_ << "\n"
    << "  pi+ in absorber  : " << simRun->npiPosIn  << "\n"
    << "  pi+ out world    : " << simRun->npiPosOut << "\n"
	<< "  Run time (s)     : " << dt                << "\n"
    << "  total steps      : " << simRun->steps     << "\n"
	<< "  debug feature ct : " << simRun->debugFeature << "\n"
    << "------------------------------\n" << G4endl;
    
}

