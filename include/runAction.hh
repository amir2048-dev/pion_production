#ifndef RUNACTION_HH
#define RUNACTION_HH
#include <string>
#include <chrono>
#include "G4UserRunAction.hh"
#include "SimConfig.hh"
#include "run.hh"

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
	std::string macroPath_;
	std::chrono::steady_clock::time_point t0_;
	std::string outDir_;
	std::string dateDir_;
	std::string date_;
	std::string timeStamp_;
	std::string runName_;
	std::string stamp_;
	
};





#endif
