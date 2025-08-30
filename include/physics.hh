#ifndef PHYSICS_HH
#define PHYSICS_HH
#include "G4VModularPhysicsList.hh"
#include "G4EmStandardPhysics.hh"
#include "QBBC.hh"
//#include "G4VHadronPhysics.hh" 

class MyphysicsList : public G4VModularPhysicsList
{
public:
	MyphysicsList();
	~MyphysicsList();
	//ConstructParticle();





};

#endif
