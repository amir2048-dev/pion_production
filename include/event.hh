#ifndef EVENT_HH
#define EVENT_HH

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include "runAction.hh"
#include "run.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "vector"
class MyEventAction : public G4UserEventAction
{
public:
	MyEventAction(const SimConfig& cfg);
	~MyEventAction();
	
	virtual void BeginOfEventAction(const G4Event*);
	virtual void EndOfEventAction(const G4Event*);
	
	void AddEvent();
    
    G4int GetEventID() {return fEventID;};
    
	
	
private:
    G4int fEventID;
    const SimConfig& cfg_;
    

};



#endif
