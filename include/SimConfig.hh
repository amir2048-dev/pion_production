#pragma once
#include "globals.hh"
#include "G4SystemOfUnits.hh"

struct SimConfig 
{
    //run name and description
    std::string runName = "testrun";
    std::string runDescription = "testing the output structure"; 

    // simulation paths
    std::string simOutDir = "../data"; // base output dir
    std::string cfgDir    = "../include"; // path to config dir
    std::string cfgFile   = "SimConfig.hh"; // config file name
    
    // world half size
    G4double worldX    = 50.0 * cm;
    G4double worldY    = 50.0 * cm;
    G4double worldZ    = 100.0 * cm;

    // absorber half size
    G4double absorberX = 0.5 * cm;
    G4double absorberY = 0.5 * cm;
    G4double absorberZ = 1.0 * cm;

    // world origin (center)
    G4double worldXOrigin = 0.0 * cm;
    G4double worldYOrigin = 0.0 * cm;
    G4double worldZOrigin = 0.0 * cm;

    // absorber origin (center) 
    G4double absorberXOrigin = 0.0 * cm;
    G4double absorberYOrigin = 0.0 * cm;
    G4double absorberZOrigin = 0.0 * cm;

    //absorber pixel size
    G4double pixelX    = 0.1 * mm;
    G4double pixelZ    = 0.1 * mm;
    
    // absorber Grid definition from pixel size and absorber size
    G4int    nAbsorberX= static_cast<G4int>(2*absorberX/pixelX);
    G4int    nAbsorberZ= static_cast<G4int>(2*absorberZ/pixelZ);

    // world pixel size (can be different from absorber pixel size)
    G4double worldPixelX = 1.0 * cm;
    G4double worldPixelZ = 1.0 * cm;

    // World Grid definition from pixel size and world size 
    G4int    nWorldX   = static_cast<G4int>(2*worldX/worldPixelX);
    G4int    nWorldZ   = static_cast<G4int>(2*worldZ/worldPixelX);
    
    // Spectra binning
    G4int    energyBins = 1000;      // 0..999 MeV (1 MeV/bin)
    G4int energymaxIndex = energyBins - 1;
    G4double energyMax  = 999.0 * MeV;

    // fluence map normalization options
    bool normalizeByArea = true;  // turn raw path length into length/area
    bool normalizePerPrimary = false; // divide by N primaries
    // Global guards
    G4double minStep    = 1.1e-6 * mm;
    G4double maxTime    = 1.0 * microsecond;

    // Run “modes”
    G4bool   runPiPlusMain = true; // main physics: detect pi+ at end detector
    G4bool   runConvStats  = true; // converter stats: gamma/e± flux, etc.
    G4bool   runDebug      = true; // ad-hoc prints/checks

    // field strength
    G4double fieldZ      = 0.1 * tesla;

    // field minStep
    G4double fieldMinStep = 0.005 * mm;

    // field solver parameters
    G4double deltaChord  = 1e-5;
    G4double minEpsilon  = 1e-5;
    G4double maxEpsilon  = 1e-3;


};
