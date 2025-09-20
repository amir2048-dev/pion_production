#ifndef RUN_HH
#define RUN_HH

#include "G4Run.hh"
#include "G4UserEventAction.hh"
#include "SimConfig.hh"


class Run : public G4Run
{
public:
	Run(const SimConfig& cfg);
	virtual ~Run(){};
//	void RecordEvent(const G4UserEventAction*) override;
	void Merge(const G4Run*) override;
	G4int ncompton;
	G4int nin;
	G4int nTransportation;
	G4int nconv;
	G4int nPhotoNuclear;
	G4int npiPosIn;
	G4int npiNegIn;
	G4int npiZerIn;
	G4int npiPosOut;
	G4int npiNegOut;
	G4int npiZerOut;
	G4int neventsOfPion;
	G4int nPhot;
	G4int pionflux[100][200] = {0};
	G4int pionfluxlarge[500][1000] = {0};
	G4int eflux[100][200] = {0};
	G4int gammaflux[100][200] = {0};
	G4int gammafluxover200[100][200] = {0};
	G4int gammacreation[100][200] = {0};
	G4int pionEnergyIn[1000] = {0};
	G4int gammaEnergy[1000] = {0};
	G4int pionEnergyOut[1000] = {0};
	G4int steps;
	G4int npiPosInElas;
private:
	const SimConfig& cfg_;
	





};
#endif
