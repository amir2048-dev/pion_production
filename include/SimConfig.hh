#pragma once
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

struct SimConfig 
{
    //run name and description
    std::string runName = "testing_angle_output";
    std::string runDescription = "testing the angle histogram output structure"; 

    // random seed control
    G4bool   useFixedSeed = false;     // if true, use fixed seed for reproducibility
    G4long   fixedSeed    = 1767432275;     // seed value when useFixedSeed = true

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

    // ----- Exit plane / angle scoring -----
    // master on/off
    G4bool   enableExitPlane       = true;     
    // exit plane position and half size     
    G4double exitPlaneZ            = 10.0 * cm;
    G4double exitPlaneHalfX        = 10.0 * cm;
    G4double exitPlaneHalfY        = 10.0 * cm;
    // angle binning using small-angle approximation
    // theta_x = atan2(px, pz), theta_y = atan2(py, pz)
    G4int    nAngleBinsThetaX      = 10;            // number of bins in x-direction
    G4int    nAngleBinsThetaY      = 10;            // number of bins in y-direction
    // auto-compute angle ranges from geometry, or override manually
    G4bool   angleAutoCompute      = true;           // if true, compute from absorber/exit-plane geometry
    // manual overrides (only used if angleAutoCompute = false)
    G4double angleMinThetaX        = -0.5;           // manual override: min angle in radians
    G4double angleMaxThetaX        = 0.5;            // manual override: max angle in radians
    G4double angleMinThetaY        = -0.5;           // manual override: min angle in radians
    G4double angleMaxThetaY        = 0.5;            // manual override: max angle in radians
    // angle selection options
    G4bool   angleIncludeBackground = false;          // if false, skip bg hist
    
    // primary gun parameters
    std::string gunParticle = "e-";

    // gun position (center) and direction
    G4ThreeVector gunPos = G4ThreeVector(0, 0, -2.0 * cm);
    G4ThreeVector gunDir = G4ThreeVector(0, 0, 1);        

    // beam spot model
    enum class BeamModel { TopHatDisk, Gaussian };
    BeamModel beamModel = BeamModel::TopHatDisk;

    // Top-hat radius (σ ignored), Gaussian σ (radius ignored) — both in transverse X/Y
    G4double beamRadius = 0.3 * cm;   // for TopHatDisk
    G4double beamSigma   = 0.3 * cm;  // for Gaussian

    // energy model
    enum class EnergyModel { Mono, Uniform, MaxwellLike };
    EnergyModel energyModel = EnergyModel::MaxwellLike;

    // energy parameters 
    G4double Emono      = 300.0 * MeV;   // Mono
    G4double EunifMin   = 200.0 * MeV;   // Uniform
    G4double Emax       = 1000.0 * MeV;
    G4double Tmaxwell   = 147.28 * MeV;  // scale in your current sampler
    G4double EminCutoff = 200.0 * MeV;   // hard low-cut

    //absorber pixel size
    G4double pixelX    = 0.1 * mm;
    G4double pixelY    = 0.1 * mm;
    G4double pixelZ    = 0.1 * mm;
        
    // absorber Grid definition from pixel size and absorber size
    G4int    nAbsorberX= static_cast<G4int>(2*absorberX/pixelX);
    G4int    nAbsorberY= static_cast<G4int>(2*absorberY/pixelY);
    G4int    nAbsorberZ= static_cast<G4int>(2*absorberZ/pixelZ);

    // world pixel size (can be different from absorber pixel size)
    G4double worldPixelX = 1.0 * cm;
    G4double worldPixelZ = 1.0 * cm;

    // World Grid definition from pixel size and world size 
    G4int    nWorldX   = static_cast<G4int>(2*worldX/worldPixelX);
    G4int    nWorldZ   = static_cast<G4int>(2*worldZ/worldPixelX);
    
    // Spectra binning (for physics outputs: pion/gamma energy)
    G4int    energyBins = 50;      // number of energy bins
    G4int    energymaxIndex = energyBins - 1;
    G4double energyMax  = Emax; // max energy for spectra histograms = primary max energy
    
    // Generator binning (hardcoded, fine resolution for generator spectrum)
    G4int    generatorEnergyBins = 1000;  // 1 MeV bins up to 1000 MeV
    
    // Helper function to get energy bin width (computed, not stored)
    G4double getEnergyBinWidth() const { return energyMax / static_cast<G4double>(energyBins); }

    // fluence map normalization options
    G4bool normalizeByArea = true;  // turn raw path length into length/area
    G4bool normalizePerPrimary = true; // divide by N primaries
    // Global guards
    G4double minStep    = 1.1e-6 * mm;
    G4double maxTime    = 1.0 * microsecond;

    // Run “modes”
    G4bool   runPiPlusMain = true; // main physics: detect pi+ at end detector
    G4bool   runConvStats  = true; // converter stats: gamma/e± flux, etc.
    G4bool   runDebug      = true; // ad-hoc prints/checks
    G4bool   runWorldMap   = true; // record world-wide pion fluence map 

    // field on/off and strength
    G4bool   enableMagneticField = false;  // if false, no magnetic field
    G4double fieldZ      = 0.1 * tesla;

    // field minStep
    G4double fieldMinStep = 0.005 * mm;

    // field solver parameters
    G4double deltaChord  = 1e-5;
    G4double minEpsilon  = 1e-5;
    G4double maxEpsilon  = 1e-3;


};
