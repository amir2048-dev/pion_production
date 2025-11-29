#include "action.hh"
#include "construction.hh"

MyActionInitialization::~MyActionInitialization()
{}
void MyActionInitialization::BuildForMaster() const
{
	MyRunAction *runAction = new MyRunAction(cfg_);
	runAction->SetMacroPath(macroPath_);
	SetUserAction(runAction);
}
void MyActionInitialization::Build() const
{
	
	MyRunAction *runAction = new MyRunAction(cfg_);
	runAction->SetMacroPath(macroPath_);
	SetUserAction(runAction);

	MyPrimaryGenerator *generator = new MyPrimaryGenerator(cfg_, runAction);
	MyEventAction *eventAction = new MyEventAction(cfg_);

	SetUserAction(generator);
	SetUserAction(eventAction);

	auto* det = dynamic_cast<const MyDetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  	const G4VPhysicalVolume* absorberPV = det ? det->GetAbsorberPV() : nullptr;
	const G4VPhysicalVolume* worldPV    = det ? det->GetWorldPV()    : nullptr;

	MySteppingAction *steppingAction = new MySteppingAction(eventAction,runAction,cfg_,absorberPV, worldPV);
	SetUserAction(steppingAction);
	
}
