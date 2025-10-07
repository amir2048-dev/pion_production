#include "event.hh"

MyEventAction::MyEventAction(const SimConfig& cfg): cfg_(cfg)
{
    fEventID = 0;
	
}
MyEventAction::~MyEventAction()
{
	
}

void MyEventAction::BeginOfEventAction(const G4Event*)
{	
	fEventID++;
}

void MyEventAction::EndOfEventAction(const G4Event*)
{
	
}


