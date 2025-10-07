#include "runAction.hh"
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <sstream>
#include <ctime>
#include <chrono>
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

	// write flattened 2D (row-major) as CSV with nx columns, nz rows
	template<typename T>
	void write_csv_2d(const std::string& path, const std::vector<T>& data, int nx, int nz, double scale = 1.0) 
	{
		if (data.empty() || nx <= 0 || nz <= 0) return;
		std::ofstream ofs(path);
		ofs.setf(std::ios::scientific, std::ios::floatfield);
		ofs << std::setprecision(3);
		for (int iz = 0; iz < nz; ++iz)
		{
			const int rowOff = iz * nx;
			for (int ix = 0; ix < nx; ++ix) 
			{
				if (ix) ofs << ',';
				ofs << (static_cast<double>(data[rowOff + ix]) * scale);
			}
			ofs << "\n";
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
	t0_ = std::chrono::steady_clock::now();

	G4cout << "[RunAction] Output dir: " << outDir_ << G4endl;
}
void MyRunAction::EndOfRunAction(const G4Run* run)
{
    if (!IsMaster()) return;
	const Run* simRun = static_cast<const Run*>(run);
	// Stop timer
	const auto t1 = std::chrono::steady_clock::now();
	const auto dt = std::chrono::duration_cast<std::chrono::seconds>(t1 - t0_).count();
    
	try_copy_file(cfg_.cfgDir + "/" + cfg_.cfgFile , outDir_ + "/" + cfg_.cfgFile); // copy config file to output dir
	try_copy_file(macroPath_, outDir_ + "/macro_used.mac"); // copy macro to output dir
	
    // normalization for fluence maps (path length / area), and optionally per-primary
  	double scaleFluence = 1.0;
  	if (cfg_.normalizeByArea)
    {
    	const double area_mm2 = (cfg_.pixelX / mm) * (cfg_.pixelZ / mm); // numeric mm^2
    	if (area_mm2 > 0.0) scaleFluence /= area_mm2;                    // length(mm)/area(mm^2) => 1/mm
	}
  	if (cfg_.normalizePerPrimary)
	{
    	const auto N = 1;//simRun->GetNumberOfEvent();// ADD FUNCTION TO GET NUMBER OF EVENTS
    	if (N > 0) scaleFluence /= static_cast<double>(N);
  	}
	if (cfg_.runPiPlusMain)
	{
		 // spectra
    	write_csv_1d(outDir_ + "/pion_energy_in.csv",  simRun->pionEnergyIn);
    	write_csv_1d(outDir_ + "/pion_energy_out.csv", simRun->pionEnergyOut);

		// absorber-local pion fluence
		write_csv_2d(outDir_ + "/pion_fluence_abs.csv",
					simRun->pionFluenceAbs, cfg_.nAbsorberX, cfg_.nAbsorberZ, scaleFluence);

		// optional world-wide map (write only if the vector is non-empty) ADD OPTION TO CFG
		if (!simRun->pionFluenceWorld.empty())
		write_csv_2d(outDir_ + "/pion_fluence_world.csv", simRun->pionFluenceWorld, cfg_.nWorldX, cfg_.nWorldZ, scaleFluence);
	}
	// Converter stats (electrons/gammas)
 	if (cfg_.runConvStats)
   	{
		// gamma spectrum inside absorber (your “gamma_energy_in”)
		write_csv_1d(outDir_ + "/gamma_energy_in.csv", simRun->gammaEnergy);

		// absorber-local fluence maps for e- and gamma
		write_csv_2d(outDir_ + "/e_fluence_abs.csv",
					simRun->eFluenceAbs, cfg_.nAbsorberX, cfg_.nAbsorberZ, scaleFluence);

		write_csv_2d(outDir_ + "/gamma_fluence_abs.csv",
					simRun->gammaFluenceAbs, cfg_.nAbsorberX, cfg_.nAbsorberZ, scaleFluence);

		write_csv_2d(outDir_ + "/gamma_fluence_over200_abs.csv",
					simRun->gammaFluenceOver200Abs, cfg_.nAbsorberX, cfg_.nAbsorberZ, scaleFluence);

		// NEED TO UPDATE TO ADD GAMMA CREATION LOCATION
		write_csv_2d(outDir_ + "/gamma_creation_abs.csv", simRun->gammaCreationAbs, cfg_.nAbsorberX, cfg_.nAbsorberZ, scaleFluence);
  	}

	std::ofstream txt(outDir_ + "/results.txt");
    txt.setf(std::ios::scientific, std::ios::floatfield);
    txt << std::setprecision(1);

    // Header
    txt << "Run folder: " << stamp_ << "\n";
    txt << "Date      : " << date_ << "\n\n";

    // From cfg
    txt << "[Config]\n";
    txt << "  run name           : " << runName_ << "\n";
    txt << "  run Description    : " << cfg_.runDescription << "\n";
    txt << "  world (half, cm)   : "
        << (cfg_.worldX/cm) << " x " << (cfg_.worldY/cm) << " x " << (cfg_.worldZ/cm) << "\n";
    txt << "  absorber(half, cm) : "
        << (cfg_.absorberX/cm) << " x " << (cfg_.absorberY/cm) << " x " << (cfg_.absorberZ/cm) << "\n";
    txt << "  field strength (T) : " << (cfg_.fieldZ/tesla) << "\n";
    // TODO: field location (add to SimConfig later)
    // txt << "  field location    : " << cfg_.fieldLocation << "\n";
    txt << "  pixel (mm)         : " << (cfg_.pixelX/mm) << " x " << (cfg_.pixelZ/mm) << "\n";
    // TODO: generator params (add to SimConfig later)
    // txt << "  generator params  : " << cfg_.generatorDesc << "\n";
    txt << "\n";

    // From macro (runtime parameters)
    txt << "[Macro]\n";
    txt << "  number of runs (N) : " << run->GetNumberOfEvent() << "\n";
    txt << "  number of threads  : " << G4RunManager::GetRunManager()->GetNumberOfThreads() << "\n";
    if (!macroPath_.empty())
      txt << "  macro file       : " << macroPath_ << "\n";
    txt << "\n";

    // From simulation counters
    txt << "[Results]\n";
    txt << "  pi+ in absorber    : " << simRun->npiPosIn << "\n";
    txt << "  pi+ outside world  : " << simRun->npiPosOut << "\n";
    txt << "\n";

    // Timing
    const auto nEvt = static_cast<double>(run->GetNumberOfEvent());
    const double perEvt = (nEvt > 0.0) ? (dt / nEvt) : 0.0;
    txt << "[Timing]\n";
    txt << "  runtime (s)        : " << dt << "\n";
    txt << "  runtime per run (s): " << perEvt << "\n";
    txt << "\n";



	G4cout
    << "\n----- End of Run (master) -----\n"
	<< "  Output dir       : " << outDir_ << "\n"
    << "  pi+ in absorber  : " << simRun->npiPosIn  << "\n"
    << "  pi+ out world    : " << simRun->npiPosOut << "\n"
	<< "  Run time (s)     : " << dt                << "\n"
    << "  total steps      : " << simRun->steps     << "\n"
    << "------------------------------\n" << G4endl;
    
}

