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
	const G4VPhysicalVolume* exitPlanePV = det ? det->GetExitPlanePV() : nullptr;
	// Get computed angle ranges from detector construction and pass to run action
	G4double angleMinThetaX = det ? det->GetAngleMinThetaX() : cfg_.angleMinThetaX;
	G4double angleMaxThetaX = det ? det->GetAngleMaxThetaX() : cfg_.angleMaxThetaX;
	G4double angleMinThetaY = det ? det->GetAngleMinThetaY() : cfg_.angleMinThetaY;
	G4double angleMaxThetaY = det ? det->GetAngleMaxThetaY() : cfg_.angleMaxThetaY;
	runAction->SetAngleRanges(angleMinThetaX, angleMaxThetaX, angleMinThetaY, angleMaxThetaY);

	MySteppingAction *steppingAction = new MySteppingAction(eventAction,runAction,cfg_,
															absorberPV, worldPV, exitPlanePV);
	SetUserAction(steppingAction);
	
}
