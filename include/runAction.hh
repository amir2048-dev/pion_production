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
	void SetAngleRanges(G4double minThetaX, G4double maxThetaX, 
	                    G4double minThetaY, G4double maxThetaY) {
		fAngleMinThetaX = minThetaX;
		fAngleMaxThetaX = maxThetaX;
		fAngleMinThetaY = minThetaY;
		fAngleMaxThetaY = maxThetaY;
	}
	G4double GetAngleMinThetaX() const { return fAngleMinThetaX; }
	G4double GetAngleMaxThetaX() const { return fAngleMaxThetaX; }
	G4double GetAngleMinThetaY() const { return fAngleMinThetaY; }
	G4double GetAngleMaxThetaY() const { return fAngleMaxThetaY; }
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
	// Computed angle ranges from detector construction
	G4double fAngleMinThetaX = 0.0;
	G4double fAngleMaxThetaX = 0.0;
	G4double fAngleMinThetaY = 0.0;
	G4double fAngleMaxThetaY = 0.0;
	
};





#endif
