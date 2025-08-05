#include "construction.hh"
#include "G4UserLimits.hh"
#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4ChordFinder.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ClassicalRK4.hh"
#include "G4Mag_UsualEqRhs.hh"

MyDetectorConstruction::MyDetectorConstruction()
{}
MyDetectorConstruction::~MyDetectorConstruction()
{}

G4VPhysicalVolume *MyDetectorConstruction::Construct()
{
	G4MagneticField* magField = new G4UniformMagField(G4ThreeVector(0.,0.,0.1*tesla));
	G4FieldManager* globalFieldMgr = G4TransportationManager::GetTransportationManager()-> GetFieldManager();
	G4Mag_UsualEqRhs* equation = new G4Mag_UsualEqRhs(magField);
	G4MagIntegratorStepper* stepper = new G4ClassicalRK4(equation);
	G4double minStep = 0.005 * mm;  // Minimum step size
	G4ChordFinder* chordFinder = new G4ChordFinder(magField, minStep, stepper);
	chordFinder->SetDeltaChord(1e-5 * mm);
	globalFieldMgr->SetChordFinder(chordFinder);
	globalFieldMgr->SetMinimumEpsilonStep(1e-5 * mm);
	globalFieldMgr->SetMaximumEpsilonStep(1e-3 * mm);
	globalFieldMgr->SetDetectorField(magField);
	G4double fact = 1;
	G4NistManager *nist = G4NistManager::Instance();
	G4Material *worldMat = nist->FindOrBuildMaterial("G4_AIR");
	//G4Box *solidWorld = new G4Box("solidWorld", 0.5*m,0.5*m,0.5*m); 
	G4Box *solidWorld = new G4Box("solidWorld", 100*cm,100*cm,200*cm*fact); // (name halfx halfy halfz)
	G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld,worldMat,"logicWorld"); // (solidworld, matirial, name)
	G4double maxStep = 0.1*mm*fact;
  	fStepLimit = new G4UserLimits(maxStep);
  	logicWorld->SetUserLimits(fStepLimit);
  	logicWorld->SetFieldManager(globalFieldMgr, true);
	G4VPhysicalVolume *physWorld = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),  logicWorld,"physWorld",0,false,0,true);
	
	G4Material *absorberMat = nist->FindOrBuildMaterial("G4_Pb"); //main code
	//G4Material *absorberMat = nist->FindOrBuildMaterial("G4_W");
	//G4Material *absorberMat = nist->FindOrBuildMaterial("G4_AIR");
	//G4Tubs 	*solidAbsorber = new G4Tubs("solidAbsorber",0., 4.3/2*cm,25./2*cm,0,2.0*pi); //name innerRad Outerrad length start angel delta angle
	
	G4Box 	*solidAbsorber = new G4Box("solidAbsorber",0.5*cm, 0.5*cm,1*cm*fact); //(name halfx halfy halfz) 
	G4LogicalVolume *logicAbsorber = new G4LogicalVolume(solidAbsorber,absorberMat,"logicAbsorber");
  	logicAbsorber->SetUserLimits(fStepLimit);
	G4VPhysicalVolume *physAbsorber = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),  logicAbsorber,"physAbsorber",logicWorld,false,0,true);
	
	//G4Box 	*solidAbsorber = new G4Box("solidAbsorber",0.1*mm, 1*cm,0.1*mm); //(name halfx halfy halfz) 
	//G4LogicalVolume *logicAbsorber = new G4LogicalVolume(solidAbsorber,absorberMat,"logicAbsorber");
	//for (G4int i=0; i<100; i++)
	//{
	//	for (G4int j=0; j<100; j++)
	//	{
	//		G4VPhysicalVolume *physAbsorber = new G4PVPlacement(0, G4ThreeVector(i*0.2*mm-1*cm,0,j*0.2*mm - 1*cm),  logicAbsorber,"physAbsorber",logicWorld,false,0,false);
	//	}
	//}
	
	return physWorld;
}

