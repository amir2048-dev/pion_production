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
	
	void AddEdep(G4double edep); 
	void AddEvent();
    
    G4int GetEventID() {return fEventID;};
    G4bool ifPionPProduced;
    G4bool ifPionNProduced;
    G4bool ifPionZProduced;
    G4bool pionfluxboolian[100][200]={0};
    G4bool pionfluxlargeboolian[500][1000]={0};
    G4bool efluxboolian[100][200]={0};
    G4bool gammafluxboolian[100][200]={0};
    G4bool gammafluxover200boolian[100][200]={0};
    G4bool gammacreationboolian[100][200]={0};
	
	
private:
    G4double fEdep;
    G4int fEventID;
    int numRows;
    int numCols;
    const SimConfig& cfg_;
    

};



#endif
