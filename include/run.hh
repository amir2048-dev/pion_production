#ifndef RUN_HH
#define RUN_HH

#include "G4Run.hh"
#include "G4UserEventAction.hh"
#include <vector>
#include "SimConfig.hh"


class Run : public G4Run
{
public:
	Run(const SimConfig& cfg);
	virtual ~Run(){};
	void Merge(const G4Run*) override;
    // scalars
  	G4int npiPosIn = 0, npiPosOut = 0, steps = 0, debugFeature = 0;
	// spectra 
  	std::vector<G4int> pionEnergyIn;   // size = cfg_.energyBins
  	std::vector<G4int> pionEnergyOut;  // size = cfg_.energyBins
  	std::vector<G4int> gammaEnergy;    // size = cfg_.energyBins
	std::vector<G4int> genratorEnergy; // size = cfg_.energyBins
	// fluence maps (accumulate path length)
	// absorber-local grids (X–Z)
	std::vector<double> pionFluenceAbs;        // length sum per cell
	std::vector<double> eFluenceAbs;
	std::vector<double> gammaFluenceAbs;
	std::vector<double> gammaFluenceOver200Abs;
	std::vector<double> gammaCreationAbs;      // if you want path length of "created" gammas' first step; or drop if not meaningful
	// world grid (X–Z)
	std::vector<double> pionFluenceWorld;
	// exit plane angle histograms
	std::vector<double> pionExitPlaneAngleHistograms;
	std::vector<double> backgroundExitPlaneAngleHistograms;


	
	// genrator beam grid (X-Y)
	std::vector<double> genratorBeamXY;

	// flatteners
  	inline int AbsIndex(int ix, int iz)   const { return ix + iz * cfg_.nAbsorberX; }
  	inline int WorldIndex(int ix, int iz) const { return ix + iz * cfg_.nWorldX;   }
	inline int BeamIndex(int ix, int iy)  const { return ix + iy * cfg_.nAbsorberX;   }

	
private:
	const SimConfig& cfg_;
	





};
#endif
