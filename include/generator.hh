#ifndef GENERATOR_HH
#define GENERATOR_HH
#include "G4VUserPrimaryGeneratorAction.hh"
//#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "SimConfig.hh"

class MyPrimaryGenerator : public G4VUserPrimaryGeneratorAction
{
public:
	MyPrimaryGenerator(const SimConfig& cfg);
	~MyPrimaryGenerator();
	virtual void GeneratePrimaries(G4Event*);

private: 
	G4ParticleGun *fParticleGun;
	const SimConfig& cfg_;
};


#endif
