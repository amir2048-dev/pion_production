#include "MyG4StepLimiterPhysics.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4StepLimiter.hh"
#include "G4UserSpecialCuts.hh"
#include "G4PhysicsConstructorFactory.hh"
G4_DECLARE_PHYSCONSTR_FACTORY(MyG4StepLimiterPhysics);
MyG4StepLimiterPhysics::MyG4StepLimiterPhysics(const G4String& name):  G4VPhysicsConstructor(name)
{;
}

MyG4StepLimiterPhysics::~MyG4StepLimiterPhysics()
{;
}
void MyG4StepLimiterPhysics::ConstructParticle() 
{}


void MyG4StepLimiterPhysics::ConstructProcess()
{
   auto aParticleIterator = GetParticleIterator();
   aParticleIterator->reset();

   G4StepLimiter* fStepLimiter = new G4StepLimiter();
   G4UserSpecialCuts* fUserSpecialCuts = new G4UserSpecialCuts();
   while ((*aParticleIterator)())
   {
     G4ParticleDefinition* particle = aParticleIterator->value();
     G4ProcessManager* pmanager = particle->GetProcessManager();
 
     if(particle->GetParticleName() == "pi+" || particle->GetParticleName() == "pi-")// || particle->GetParticleName() == "gamma" || particle->GetParticleName() == "e-") 
     {
     	pmanager->AddDiscreteProcess(fStepLimiter);
     	pmanager->AddDiscreteProcess(fUserSpecialCuts);
     }
     else if(particle->GetParticleName() == "e-" || particle->GetParticleName() == "e+" || particle->GetParticleName() == "gamma")
     {
     	pmanager->AddDiscreteProcess(fStepLimiter); 
     	pmanager->AddDiscreteProcess(fUserSpecialCuts);
     }
     
   }
}
 
