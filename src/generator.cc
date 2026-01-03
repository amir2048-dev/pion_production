#include "generator.hh"
#include "Randomize.hh"
#include "G4ParticleDefinition.hh"
#include "G4ios.hh"
#include <cmath>

inline double f_shape(double E, double T) 
{
    if (E <= 0.0 || T <= 0.0) return 0.0;
    return std::sqrt(E) * std::pow(T, -1.5) * std::exp(-E / T);
}

MyPrimaryGenerator::MyPrimaryGenerator(const SimConfig& cfg, MyRunAction* runAction ): cfg_(cfg), runAction_(runAction)
{
	particleGun_ = new G4ParticleGun(1);
	
	// prticle
	G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
	auto* def   = particleTable->FindParticle(cfg_.gunParticle);
	if (!def)
	{
    	G4Exception("MyPrimaryGenerator", "GEN001", JustWarning, ("Unknown particle name: " + cfg_.gunParticle).c_str());
    	def = particleTable->FindParticle("e-");
    }
  	particleGun_->SetParticleDefinition(def);

	 // direction (normalized; if zero, fall back to +z)
  	G4ThreeVector dir = cfg_.gunDir;
  	if (dir.mag2() == 0) dir = G4ThreeVector(0,0,1);
  	particleGun_->SetParticleMomentumDirection(dir.unit());
	particleGun_->SetParticlePosition(cfg_.gunPos);
}

MyPrimaryGenerator::~MyPrimaryGenerator()
{
	delete particleGun_; 
}

void MyPrimaryGenerator::GeneratePrimaries(G4Event *anEvent)
{
	// Sample position in X/Y at gun Z
	const auto xy = samplePositionXY_();
  	const auto gunPos = G4ThreeVector(cfg_.gunPos.x() + xy.x(), cfg_.gunPos.y() + xy.y(), cfg_.gunPos.z());
  	particleGun_->SetParticlePosition(gunPos);

	// increment generator beam histogram
	const auto x = gunPos.x();
	const auto y = gunPos.y();

	// map to indices with world origin at (x0,y0)
	const double gx = (x - cfg_.gunPos.x()) / cfg_.pixelX + 0.5 * cfg_.nAbsorberX;
	const double gy = (y - cfg_.gunPos.y()) / cfg_.pixelY + 0.5 * cfg_.nAbsorberY;
	int ix = static_cast<int>(std::floor(gx));
	int iy = static_cast<int>(std::floor(gy));

	if (0 <= ix && ix < cfg_.nAbsorberX && 0 <= iy && iy < cfg_.nAbsorberY)
	{
		const int k = runAction_->fRun->BeamIndex(ix, iy);
		runAction_->fRun->genratorBeamXY[k] += 1.0;   // count primaries per pixel
	}
	
	// Sample energy (apply cutoff + model)
  	G4double E = sampleEnergy_();
  	particleGun_->SetParticleEnergy(E);
	// increment generator energy histogram (1 MeV bins)
	int bin = static_cast<int>(E / MeV);
	if (bin < 0) bin = 0;
	if (bin >= cfg_.generatorEnergyBins) bin = cfg_.generatorEnergyBins - 1;
	runAction_->fRun->genratorEnergy[bin]++;

  	particleGun_->GeneratePrimaryVertex(anEvent);
}
G4ThreeVector MyPrimaryGenerator::samplePositionXY_() const
{
	switch (cfg_.beamModel)
  	{
		case SimConfig::BeamModel::TopHatDisk:
		{
	  		// uniform in circle of radius cfg_.beamRadius
	  		// Uniform disk: r = R * sqrt(u), theta = 2πv
			const G4double R = cfg_.beamRadius;
			const G4double u = G4UniformRand();
			const G4double v = G4UniformRand();
			const G4double r = R * std::sqrt(u);
			const G4double th = 2.0 * CLHEP::pi * v;
			return { r * std::cos(th), r * std::sin(th), 0.0 };
		}
		case SimConfig::BeamModel::Gaussian:
		{
	  		// 2D Gaussian in X/Y with sigma=cfg_.beamSigma
	  		G4double x = G4RandGauss::shoot(0, cfg_.beamSigma);
	  		G4double y = G4RandGauss::shoot(0, cfg_.beamSigma);
	  		return G4ThreeVector(x, y, 0);
		}
		default:
	  	return G4ThreeVector(0, 0, 0); // should not happen
  	}
}
G4double MyPrimaryGenerator::sampleEnergy_() const 
{
	using EM = SimConfig::EnergyModel;
	switch (cfg_.energyModel)
	{
		case EM::Mono:
		return cfg_.Emono;

		case EM::Uniform: 
		{
			const G4double a = cfg_.EunifMin;
			const G4double b = cfg_.Emax;
			const G4double E = a + (b - a) * G4UniformRand();
			return (E < cfg_.EminCutoff) ? cfg_.EminCutoff : E;
		}

		case EM::MaxwellLike:
		default:
		return sampleMaxwellLike_(cfg_.Tmaxwell, cfg_.EminCutoff, cfg_.Emax);
	}
}

G4double MyPrimaryGenerator::sampleMaxwellLike_(G4double T, G4double Emin, G4double Emax)
{
    if (T <= 0.0) return std::max(Emin, 1.0*MeV);
    if (Emax <= Emin) return Emin;

    const double Epeak = 0.5 * T;
    const double Em    = std::min(std::max(Emin, Epeak), Emax); 
    const double ymax  = f_shape(Em, T);
    if (ymax <= 0.0)   return Emin;

    while (true)
	{
        const double E = Emin + (Emax - Emin) * G4UniformRand();
        const double y = ymax * G4UniformRand();
        if (y <= f_shape(E, T)) return E;
    }
}

	