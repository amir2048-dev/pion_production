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

MyPrimaryGenerator::MyPrimaryGenerator(const SimConfig& cfg): cfg_(cfg)
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
	const auto xy = samplePositionXY_();
  	const auto gunPos = G4ThreeVector(cfg_.gunPos.x() + xy.x(), cfg_.gunPos.y() + xy.y(), cfg_.gunPos.z());
  	particleGun_->SetParticlePosition(gunPos);

	// Sample energy (apply cutoff + model + scale)
  	G4double E = sampleEnergy_();
  	E *= cfg_.energyScale;
  	particleGun_->SetParticleEnergy(E);
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

	