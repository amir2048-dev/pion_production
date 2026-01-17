#include "ElectromagneticPhysics.hh"

#include "G4BuilderType.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicsListHelper.hh"

#include "G4EmParameters.hh"
#include "G4LossTableManager.hh"
#include "G4UAtomicDeexcitation.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4RayleighScattering.hh"

#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuMultipleScattering.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"

#include "G4ionIonisation.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4NuclearStopping.hh"

#include "G4SystemOfUnits.hh"

ElectromagneticPhysics::ElectromagneticPhysics(const G4String& name)
: G4VPhysicsConstructor(name)
{
  SetPhysicsType(bElectromagnetic);

  auto* param = G4EmParameters::Instance();
  param->SetDefaults();

  // ------------------------------
  // DEV (fast iteration) settings
  // ------------------------------
  // Looser step functions -> fewer steps -> big speedup in thick Pb showers.
  // These primarily affect tracking granularity (not which processes exist).
  param->SetStepFunction(0.3, 1*mm);             // e±  (DEV)
  param->SetStepFunctionMuHad(0.3, 1*mm);        // mu/hadrons (DEV)
  param->SetStepFunctionLightIons(0.3, 1*mm);    // light ions (DEV)
  param->SetStepFunctionIons(0.3, 1*mm);         // ions (DEV)

  // Disable atomic de-excitation in DEV (saves low-energy secondaries)
  // Precision impact: removes fluorescence/Auger X-rays and very low-E electrons.
  param->SetFluo(false);
  param->SetAuger(false);
  param->SetPixe(false);

  // If you later enable de-excitation, this controls whether cuts are ignored.
  // Keeping IgnoreCut=true can increase low-energy production (slower).
  param->SetDeexcitationIgnoreCut(false);        // DEV: respect cuts (or just keep fluo off)

  // ------------------------------
  // FINAL (more precise) options
  // ------------------------------
  // Tighten step functions for final detector/noise studies:
  // param->SetStepFunction(0.2, 100*um);         // your original e±
  // param->SetStepFunctionMuHad(0.1, 10*um);
  // param->SetStepFunctionLightIons(0.1, 10*um);
  // param->SetStepFunctionIons(0.1, 1*um);
  //
  // Enable atomic de-excitation when you care about low-energy X-rays/Auger:
  // param->SetFluo(true);
  // param->SetAuger(true);
  // param->SetPixe(true);
  // param->SetDeexcitationIgnoreCut(true);
}

ElectromagneticPhysics::~ElectromagneticPhysics() {}

void ElectromagneticPhysics::ConstructProcess()
{
  auto* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  auto particleIterator = GetParticleIterator();
  particleIterator->reset();

  while ((*particleIterator)()) {
    G4ParticleDefinition* particle = particleIterator->value();
    const G4String& particleName = particle->GetParticleName();

    if (particleName == "gamma") {

      // DEV: Rayleigh is mostly low-energy; you can disable for speed.
      // Precision impact: small at MeV–GeV, matters at keV–100s keV.
      // ph->RegisterProcess(new G4RayleighScattering,  particle);  // FINAL option

      ph->RegisterProcess(new G4PhotoElectricEffect, particle);
      ph->RegisterProcess(new G4ComptonScattering,   particle);
      ph->RegisterProcess(new G4GammaConversion,     particle);

    } else if (particleName == "e-") {

      ph->RegisterProcess(new G4eMultipleScattering(), particle);
      ph->RegisterProcess(new G4eIonisation,           particle);
      ph->RegisterProcess(new G4eBremsstrahlung(),     particle);

    } else if (particleName == "e+") {

      ph->RegisterProcess(new G4eMultipleScattering(), particle);
      ph->RegisterProcess(new G4eIonisation,           particle);
      ph->RegisterProcess(new G4eBremsstrahlung(),     particle);
      ph->RegisterProcess(new G4eplusAnnihilation(),   particle);

    } else if (particleName == "mu+" || particleName == "mu-") {

      ph->RegisterProcess(new G4MuMultipleScattering(), particle);
      ph->RegisterProcess(new G4MuIonisation,           particle);
      ph->RegisterProcess(new G4MuBremsstrahlung(),     particle);
      ph->RegisterProcess(new G4MuPairProduction(),     particle);

    } else if (particleName == "proton" || particleName == "pi-" || particleName == "pi+") {

      ph->RegisterProcess(new G4hMultipleScattering(), particle);
      ph->RegisterProcess(new G4hIonisation,           particle);

    } else if (particleName == "alpha" || particleName == "He3") {

      ph->RegisterProcess(new G4hMultipleScattering(), particle);
      ph->RegisterProcess(new G4ionIonisation,         particle);
      ph->RegisterProcess(new G4NuclearStopping(),     particle);

    } else if (particleName == "GenericIon") {

      ph->RegisterProcess(new G4hMultipleScattering(), particle);
      auto* ionIoni = new G4ionIonisation();
      ionIoni->SetEmModel(new G4IonParametrisedLossModel());
      ph->RegisterProcess(ionIoni,                  particle);
      ph->RegisterProcess(new G4NuclearStopping(),  particle);

    } else if ((!particle->IsShortLived()) &&
               (particle->GetPDGCharge() != 0.0) &&
               (particleName != "chargedgeantino")) {

      ph->RegisterProcess(new G4hMultipleScattering(), particle);
      ph->RegisterProcess(new G4hIonisation(),         particle);
    }
  }

  // Atomic de-excitation object:
  // DEV: if Fluo/Auger/PIXE are all false, installing this is not needed.
  // FINAL: enable if you turn on fluo/auger/pixe.
  //
  // auto* de = new G4UAtomicDeexcitation();
  // G4LossTableManager::Instance()->SetAtomDeexcitation(de);
}
