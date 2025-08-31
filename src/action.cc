#include "action.hh"
#include "construction.hh"

MyActionInitialization::MyActionInitialization()
{}
MyActionInitialization::~MyActionInitialization()
{}
void MyActionInitialization::BuildForMaster() const
{
	MyRunAction *runAction = new MyRunAction();
	SetUserAction(runAction);
}
void MyActionInitialization::Build() const
{
	MyPrimaryGenerator *generator = new MyPrimaryGenerator();
	MyRunAction *runAction = new MyRunAction();
	MyEventAction *eventAction = new MyEventAction();

	SetUserAction(generator);
	SetUserAction(runAction);
	SetUserAction(eventAction);

	auto* det = dynamic_cast<const MyDetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  	const G4VPhysicalVolume* absorberPV = det ? det->GetAbsorberPV() : nullptr;
	const G4VPhysicalVolume* worldPV    = det ? det->GetWorldPV()    : nullptr;

	MySteppingAction *steppingAction = new MySteppingAction(eventAction,runAction,absorberPV, worldPV);
	SetUserAction(steppingAction);
	
}
