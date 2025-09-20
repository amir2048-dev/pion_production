#ifndef RUNACTION_HH
#define RUNACTION_HH
#include "G4UserRunAction.hh"
#include "G4AnalysisManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Run.hh"
#include "run.hh"
#include "SimConfig.hh"
#include <mutex>

class MyRunAction : public G4UserRunAction
{
public: 
	MyRunAction(const SimConfig& cfg): cfg_(cfg) {};
	~MyRunAction();
	virtual G4Run* GenerateRun();
	virtual void BeginOfRunAction(const G4Run*);
	virtual void EndOfRunAction(const G4Run*);
	Run* fRun;
	
private:
	const SimConfig& cfg_;
	
};





#endif
