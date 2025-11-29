#ifndef GENERATOR_HH
#define GENERATOR_HH
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4ThreeVector.hh"
#include "SimConfig.hh"
#include "runAction.hh"

class MyPrimaryGenerator : public G4VUserPrimaryGeneratorAction
{
public:
	MyPrimaryGenerator(const SimConfig& cfg, MyRunAction* runAction);
	~MyPrimaryGenerator();
	virtual void GeneratePrimaries(G4Event*);

private: 
	G4ParticleGun* particleGun_ = nullptr;
	MyRunAction* runAction_ = nullptr;
	const SimConfig& cfg_;
	// samplers
  	G4ThreeVector samplePositionXY_() const;   // transverse spot at gun Z
  	G4double      sampleEnergy_() const;       // energy in [MeV] with model+cutoff

  	// helpers
  	static G4double sampleMaxwellLike_(G4double T, G4double Emin, G4double Emax);
};


#endif
