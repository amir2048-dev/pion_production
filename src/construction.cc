#include "construction.hh"
#include "G4SystemOfUnits.hh"
#include "G4AutoDelete.hh"
#include "G4UserLimits.hh"
#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4ChordFinder.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ClassicalRK4.hh"
#include "G4Mag_UsualEqRhs.hh"



MyDetectorConstruction::~MyDetectorConstruction()
{}
G4ThreadLocal G4MagneticField*         MyDetectorConstruction::tMagField  = nullptr;
G4ThreadLocal G4Mag_UsualEqRhs*        MyDetectorConstruction::tEquation  = nullptr;
G4ThreadLocal G4MagIntegratorStepper*  MyDetectorConstruction::tStepper   = nullptr;
G4ThreadLocal G4ChordFinder*           MyDetectorConstruction::tChord     = nullptr;
void MyDetectorConstruction::ConstructSDandField()
{
	if (!tMagField)
    {	
		// ---- your field ----
		tMagField = new G4UniformMagField(G4ThreeVector(0., 0., cfg_.fieldZ));
		G4AutoDelete::Register(tMagField);

		// ---- your equation & stepper (per-thread) ----
		tEquation = new G4Mag_UsualEqRhs(tMagField);
		G4AutoDelete::Register(tEquation);

		tStepper  = new G4ClassicalRK4(tEquation);
		G4AutoDelete::Register(tStepper);

		const G4double minStep = cfg_.minStep ;   
		tChord = new G4ChordFinder(tMagField, minStep, tStepper);
		G4AutoDelete::Register(tChord);

		// ---- attach to the thread-local FieldManager and tune tolerances ----
		auto fm = G4TransportationManager::GetTransportationManager()->GetFieldManager();
		fm->SetDetectorField(tMagField);
		fm->SetChordFinder(tChord);

		tChord->SetDeltaChord(cfg_.deltaChord);
		fm->SetMinimumEpsilonStep(cfg_.minEpsilon);
		fm->SetMaximumEpsilonStep(cfg_.maxEpsilon);
	}
    
}
G4VPhysicalVolume *MyDetectorConstruction::Construct()
{
	G4FieldManager* globalFieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
	G4NistManager *nist = G4NistManager::Instance();
	G4Material *worldMat = nist->FindOrBuildMaterial("G4_AIR"); 
	G4Box *solidWorld = new G4Box("solidWorld", cfg_.worldX, cfg_.worldY, cfg_.worldZ ); // (name halfx halfy halfz)
	fSolidWorld = solidWorld;
	G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld,worldMat,"logicWorld"); // (solidworld, matirial, name)
	fPhysWorld = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),  logicWorld,"physWorld",0,false,0,true);
	G4Material *absorberMat = nist->FindOrBuildMaterial("G4_Pb"); 
	G4Box 	*solidAbsorber = new G4Box("solidAbsorber",cfg_.absorberX, cfg_.absorberY, cfg_.absorberZ); //(name halfx halfy halfz) 
	fSolidAbsorber = solidAbsorber;
	G4LogicalVolume *logicAbsorber = new G4LogicalVolume(solidAbsorber,absorberMat,"logicAbsorber");
	fPhysAbsorber = new G4PVPlacement(0, G4ThreeVector(cfg_.absorberXOrigin,cfg_.absorberYOrigin,cfg_.absorberZOrigin),  logicAbsorber,"physAbsorber",logicWorld,false,0,true);
	
	return fPhysWorld;
}

