#include "physics.hh"

MyphysicsList::MyphysicsList()
{
	RegisterPhysics (new G4EmStandardPhysics());
	//RegisterPhysics (new G4VHadronPhysics());
	
}
MyphysicsList::~MyphysicsList()
{}
