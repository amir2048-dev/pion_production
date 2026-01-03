#ifndef STEPPING_HH
#define STEPPING_HH
#include "G4UserSteppingAction.hh"
#include "G4Step.hh"
#include "G4VProcess.hh"
#include "G4SystemOfUnits.hh"
#include "construction.hh"
#include "event.hh"
#include "runAction.hh"
#include "SimConfig.hh"
class MySteppingAction: public G4UserSteppingAction
{
public:
	MySteppingAction(MyEventAction* eventAction, MyRunAction* runAction,const SimConfig& cfg, 
					 const G4VPhysicalVolume* absorberPV,const G4VPhysicalVolume* worldPV, 
					 const G4VPhysicalVolume* exitPlanePV):
	fEventAction(eventAction), fRunAction(runAction), cfg_(cfg),
	fAbsorberPV(absorberPV), fWorldPV(worldPV), fExitPlanePV(exitPlanePV) {};
	~MySteppingAction();
	void UserSteppingAction(const G4Step* step) override;
	
private:
	MyEventAction* fEventAction = nullptr;
  	MyRunAction*   fRunAction   = nullptr;
  	const G4VPhysicalVolume* fAbsorberPV = nullptr; 
  	const G4VPhysicalVolume* fWorldPV    = nullptr;
	const G4VPhysicalVolume* fExitPlanePV = nullptr;
	G4bool afterPhotonuclear = false;
	const SimConfig& cfg_;

};
#endif
