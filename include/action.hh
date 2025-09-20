#ifndef ACTION_HH
#define ACTION_HH

#include "G4VUserActionInitialization.hh"
#include "generator.hh"
#include "stepping.hh"
#include "event.hh"
#include "runAction.hh" 
#include "SimConfig.hh"

class MyActionInitialization : public G4VUserActionInitialization
{
public: 
	MyActionInitialization(const SimConfig& cfg): cfg_(cfg) {};
	~MyActionInitialization();
	void Build() const override;
	void BuildForMaster() const override;
private:
	const SimConfig& cfg_;
	
};




#endif
