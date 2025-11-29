#ifndef RUNACTION_HH
#define RUNACTION_HH
#include <string>
#include <chrono>
#include "G4UserRunAction.hh"
#include "SimConfig.hh"
#include "run.hh"
#include "G4Timer.hh"

class MyRunAction : public G4UserRunAction
{
public: 
	MyRunAction(const SimConfig& cfg): cfg_(cfg) {};
	~MyRunAction();
	virtual G4Run* GenerateRun();
	virtual void BeginOfRunAction(const G4Run*);
	virtual void EndOfRunAction(const G4Run*);
	Run* fRun;
	void SetMacroPath(const std::string& path) { macroPath_ = path; }
private:
	const SimConfig& cfg_;
	G4double spectrumNorm_ = 1.0;
	std::string macroPath_;
	G4Timer myTimer_;
	std::string outDir_;
	std::string dateDir_;
	std::string date_;
	std::string timeStamp_;
	std::string runName_;
	std::string stamp_;
	
};





#endif
