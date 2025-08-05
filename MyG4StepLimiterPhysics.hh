#ifndef MyG4StepLimiterPhysics_h
#define MyG4StepLimiterPhysics_h 1
 
 
#include "G4VPhysicsConstructor.hh"
#include "globals.hh"
 
class G4StepLimiter;
class G4UserSpecialCuts;
 
 
class MyG4StepLimiterPhysics : public G4VPhysicsConstructor
{
public:
 
	MyG4StepLimiterPhysics(const G4String& name = "stepLimiter");
	virtual ~MyG4StepLimiterPhysics();
	virtual void ConstructParticle();
	virtual void ConstructProcess();
private:
	MyG4StepLimiterPhysics & operator=(const MyG4StepLimiterPhysics &right);
	MyG4StepLimiterPhysics(const MyG4StepLimiterPhysics&);
};
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 #endif
