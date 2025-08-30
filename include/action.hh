#ifndef ACTION_HH
#define ACTION_HH

#include "G4VUserActionInitialization.hh"
#include "generator.hh"
#include "stepping.hh"
#include "event.hh"
#include "runAction.hh" 

class MyActionInitialization : public G4VUserActionInitialization
{
public: 
	MyActionInitialization();
	~MyActionInitialization();
	void Build() const override;
	void BuildForMaster() const override;
	

	
};




#endif
