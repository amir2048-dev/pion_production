#ifndef CONSTRUCTION_HH
#define CONSTRUCTION_HH
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"                 // defines G4ThreadLocal, units, etc.
#include "G4MagneticField.hh"         // base magnetic field
#include "G4UniformMagField.hh"       // if you use uniform fields
#include "G4Mag_UsualEqRhs.hh"        // equation of motion
#include "G4MagIntegratorStepper.hh"  // abstract stepper
#include "G4ChordFinder.hh"           // chord finder
#include "SimConfig.hh"

class MyDetectorConstruction : public G4VUserDetectorConstruction
{
	public:
		explicit MyDetectorConstruction(const SimConfig& cfg) : cfg_(cfg) {};
		~MyDetectorConstruction();
		virtual G4VPhysicalVolume *Construct();
		virtual void ConstructSDandField() override;
		G4VPhysicalVolume* GetWorldPV()    const { return fPhysWorld; }
		G4VPhysicalVolume* GetAbsorberPV() const { return fPhysAbsorber; }
		G4VPhysicalVolume* GetExitPlanePV() const { return fPhysExitPlane; }
		const G4Box* GetWorldBox() const { return fSolidWorld; }
		const G4Box* GetAbsorberBox() const { return fSolidAbsorber; }
		const G4Box* GetExitPlaneBox() const { return fSolidExitPlane; }
		// Get computed or manual angle ranges
		G4double GetAngleMinThetaX() const { return fAngleMinThetaX; }
		G4double GetAngleMaxThetaX() const { return fAngleMaxThetaX; }
		G4double GetAngleMinThetaY() const { return fAngleMinThetaY; }
		G4double GetAngleMaxThetaY() const { return fAngleMaxThetaY; }	
	private:
		const SimConfig& cfg_;
		G4UserLimits* fStepLimit = nullptr; // pointer to user step limits
		static G4ThreadLocal G4MagneticField*     tMagField;
		static G4ThreadLocal G4Mag_UsualEqRhs*    tEquation;
		static G4ThreadLocal G4MagIntegratorStepper* tStepper;
		static G4ThreadLocal G4ChordFinder*       tChord;
		G4VPhysicalVolume* fPhysWorld = nullptr;
		G4VPhysicalVolume* fPhysAbsorber = nullptr;
		G4VPhysicalVolume* fPhysExitPlane = nullptr;
	G4Box*             fSolidWorld   = nullptr;
	G4Box*             fSolidAbsorber= nullptr;
	G4Box*             fSolidExitPlane= nullptr;
	// Computed angle ranges (set during Construct())
	G4double           fAngleMinThetaX = 0.0;
	G4double           fAngleMaxThetaX = 0.0;
	G4double           fAngleMinThetaY = 0.0;
	G4double           fAngleMaxThetaY = 0.0;

};

#endif
