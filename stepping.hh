#ifndef STEPPING_HH
#define STEPPING_HH
#include "G4UserSteppingAction.hh"
#include "G4Step.hh"
#include "G4VProcess.hh"
#include "G4SystemOfUnits.hh"
#include "construction.hh"
#include "event.hh"
#include "runAction.hh"
class MySteppingAction: public G4UserSteppingAction
{
public:
	MySteppingAction(MyEventAction* eventAction, MyRunAction* runAction);
	~MySteppingAction();
	virtual void UserSteppingAction(const G4Step*);
	
private:
	MyEventAction *fEventAction;
	MyRunAction *fRunAction;
	G4bool afterPhotonuclear;



};
#endif
