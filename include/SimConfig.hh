#pragma once
#include "globals.hh"
#include "G4SystemOfUnits.hh"

struct SimConfig 
{
    // world size
    G4double worldX    = 100.0 * cm;
    G4double worldY    = 100.0 * cm;
    G4double worldZ    = 200.0 * cm;

    // absorber size
    G4double absorberX = 1.0 * cm;
    G4double absorberY = 1.0 * cm;
    G4double absorberZ = 2.0 * cm;

    // world origin (center)
    G4double worldXOrigin = 0.0 * cm;
    G4double worldYOrigin = 0.0 * cm;
    G4double worldZOrigin = 0.0 * cm;

    // absorber origin (center) 
    G4double absorberXOrigin = 0.0 * cm;
    G4double absorberYOrigin = 0.0 * cm;
    G4double absorberZOrigin = 0.0 * cm;

    // pixel size
    G4double pixelX    = 0.1 * mm;
    G4double pixelZ    = 0.1 * mm;

    // World Grid definition from pixel size and world size 
    G4int    nWorldX   = static_cast<G4int>(worldX/pixelX);
    G4int    nWorldZ   = static_cast<G4int>(worldZ/pixelZ);
    
    // absorber Grid definition from pixel size and absorber size
    G4int    nAbsorberX= static_cast<G4int>(absorberX/pixelX);
    G4int    nAbsorberZ= static_cast<G4int>(absorberZ/pixelZ);

    // Spectra binning
    G4int    energyBins = 1000;      // 0..999 MeV (1 MeV/bin)
    G4double energyMax  = 999.0 * MeV;

    // Global guards
    G4double minStep    = 1.1e-6 * mm;
    G4double maxStep    = 0.1 * mm;
    G4double maxTime    = 1.0 * microsecond;

    // Run “modes”
    G4bool   runPiPlusMain = true;
    G4bool   runConvStats  = false;
    G4bool   runDebug      = false;

    // field strength
    G4double fieldZ      = 0.1 * tesla;

    // field minStep
    G4double fieldMinStep = 0.005 * mm;

    // field solver parameters
    G4double deltaChord  = 1e-5;
    G4double minEpsilon  = 1e-5;
    G4double maxEpsilon  = 1e-3;


};
