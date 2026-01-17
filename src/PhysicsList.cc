// PhysicsList.cc
// Clean, consistent, "light for development" baseline with clearly marked precision options.
// Goal: 300 MeV-ish e- on thick Pb -> bremsstrahlung -> photonuclear pions -> transport in B-field
// Keep flexibility for future detector/noise realism without paying the cost now.

#include "PhysicsList.hh"

#include "G4SystemOfUnits.hh"
#include "G4NuclideTable.hh"

// Hadronic / ions
#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsXS.hh"      // optional (slower, more detailed)
#include "HadronElasticPhysicsHP.hh"        // your optional HP neutron elastic (<20 MeV)

#include "G4HadronPhysicsFTFP_BERT.hh"
#include "G4HadronPhysicsFTFP_BERT_HP.hh"   // optional (slower) full HP neutrons (incl capture/inelastic)
#include "G4IonElasticPhysics.hh"
#include "G4IonPhysicsXS.hh"
#include "G4StoppingPhysics.hh"

// Gamma-nuclear / photonuclear
#include "GammaNuclearPhysics.hh"

// EM / decay
#include "ElectromagneticPhysics.hh"        // your EM physics (brems, MSC, etc.)
// #include "G4EmStandardPhysics_option3.hh" // optional: Geant4 reference EM list for validation
#include "G4DecayPhysics.hh"

// Radioactive decay (usually OFF for dev speed)
#include "RadioactiveDecayPhysics.hh"
// #include "G4RadioactiveDecayPhysics.hh"  // optional: standard Geant4 RDM

// Step limiter
#include "MyG4StepLimiterPhysics.hh"

PhysicsList::PhysicsList() : G4VModularPhysicsList()
{
  const G4int verb = 1;
  SetVerboseLevel(verb);

  // Mandatory for G4NuclideTable
  const G4double meanLife = 1 * nanosecond;
  const G4double halfLife = meanLife * std::log(2.0);
  G4NuclideTable::GetInstance()->SetThresholdOfHalfLife(halfLife);

  // ------------------------------------------------------------
  // HADRON ELASTIC (choose ONE)
  // ------------------------------------------------------------
  // DEV (recommended): standard elastic (usually fastest)
  RegisterPhysics(new G4HadronElasticPhysics(verb));

  // OPTION A (more detail, can be slower): elastic with extra XS handling
  // RegisterPhysics(new G4HadronElasticPhysicsXS(verb));

  // OPTION B (neutron elastic < 20 MeV more accurate, slower):
  // Use when neutron transport in shielding/detector background matters.
  // RegisterPhysics(new HadronElasticPhysicsHP(verb));
  //   If you enable thermal scattering in that class, it must be done pre-init:
  //   /testhadr/phys/thermalScattering true

  // ------------------------------------------------------------
  // HADRON INELASTIC (choose ONE)
  // ------------------------------------------------------------
  // DEV (recommended): coherent, standard, good performance
  RegisterPhysics(new G4HadronPhysicsFTFP_BERT(verb));

  // OPTION (slower): full HP neutron treatment (elastic+inelastic+capture etc. for low-E n)
  // Use only for final detector/noise studies where neutron backgrounds matter.
  // RegisterPhysics(new G4HadronPhysicsFTFP_BERT_HP(verb));

  // IMPORTANT: Do NOT stack other full hadron constructors (e.g., QGSP_BIC_* on top of FTFP_BERT).
  // It can introduce overlapping models/processes and harm both speed and interpretability.

  // ------------------------------------------------------------
  // IONS / STOPPING
  // ------------------------------------------------------------
  RegisterPhysics(new G4IonElasticPhysics(verb));
  RegisterPhysics(new G4IonPhysicsXS(verb));
  RegisterPhysics(new G4StoppingPhysics(verb));

  // ------------------------------------------------------------
  // GAMMA-NUCLEAR / PHOTONUCLEAR (must-have for photo-pion production)
  // ------------------------------------------------------------
  RegisterPhysics(new GammaNuclearPhysics("gamma"));

  // ------------------------------------------------------------
  // EM (must-have for bremsstrahlung and EM shower development)
  // ------------------------------------------------------------
  RegisterPhysics(new ElectromagneticPhysics());
  // OPTION (for cross-checking physics/results): Geant4 reference EM list
  // RegisterPhysics(new G4EmStandardPhysics_option3());

  // ------------------------------------------------------------
  // DECAYS
  // ------------------------------------------------------------
  RegisterPhysics(new G4DecayPhysics());

  // Radioactive decay: typically OFF for dev speed unless studying activation backgrounds
  RegisterPhysics(new RadioactiveDecayPhysics());
  // OPTION: standard Geant4 radioactive decay
  // RegisterPhysics(new G4RadioactiveDecayPhysics());

  // ------------------------------------------------------------
  // STEP LIMITER (keep if your workflow relies on it)
  // ------------------------------------------------------------
  //RegisterPhysics(new MyG4StepLimiterPhysics());
}

PhysicsList::~PhysicsList() {}

void PhysicsList::SetCuts()
{
  // Production cuts table should span a wide energy range (sanity/stability).
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(250 * eV, 100 * GeV);

  // ------------------------------------------------------------
  // DEV CUTS (fast iteration)
  // ------------------------------------------------------------
  // These are big cuts -> fewer low-energy e±/γ secondaries tracked -> much faster.
  // This can under-estimate detector noise later, so treat as "development mode".
  SetCutValue(100 * mm, "e-");
  SetCutValue(100 * mm, "e+");
  SetCutValue(50  * mm, "gamma");

  // OPTION (more realistic EM shower / detector noise): reduce global cuts
  // SetCutValue(1 * mm, "e-");
  // SetCutValue(1 * mm, "e+");
  // SetCutValue(1 * mm, "gamma");

  // BEST PRACTICE FOR FUTURE DETECTOR WORK:
  // Keep large GLOBAL cuts for speed, and apply SMALL cuts only in detector/shielding regions
  // via G4Region + G4ProductionCuts (regional cuts). That gives realistic PMT noise without
  // slowing the whole simulation.

  // If you later want to suppress very soft hadrons globally, set a nonzero proton cut
  // (leaving it unset uses default; setting 0 mm may increase secondaries/CPU).
  // SetCutValue(?? * mm, "proton");
}
