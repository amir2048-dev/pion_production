#include "run.hh"
#include "event.hh"
Run::Run(const SimConfig& cfg) : cfg_(cfg)
{
    //spectra assignment
	pionEnergyIn.assign(cfg_.energyBins, 0);
	pionEnergyOut.assign(cfg_.energyBins, 0);
	gammaEnergy.assign(cfg_.energyBins, 0);
	//flunce absorber maps assignment
	const int nA = cfg_.nAbsorberX * cfg_.nAbsorberZ;
	pionFluenceAbs.assign(nA, 0.0);
	eFluenceAbs.assign(nA, 0.0);
	gammaFluenceAbs.assign(nA, 0.0);
	gammaFluenceOver200Abs.assign(nA, 0.0);
	gammaCreationAbs.assign(nA, 0.0);
	//flunce world map assignment
	const int nW = cfg_.nWorldX * cfg_.nWorldZ;
	pionFluenceWorld.assign(nW, 0.0);		
}
void Run::Merge(const G4Run* aRun)
{
	
	const Run* localRun = static_cast<const Run*>(aRun);
	// mergeing all counters
	npiPosIn+=localRun->npiPosIn;
	npiPosOut+=localRun->npiPosOut;
	steps+=localRun->steps;
	// mergeing pion/e/gamma location in the absorber
	const G4int nA = cfg_.nAbsorberX * cfg_.nAbsorberZ;
	for (int i=0;i<nA;i++)
	{
		pionFluenceAbs[i]+=localRun->pionFluenceAbs[i];
		eFluenceAbs[i]+=localRun->eFluenceAbs[i];
		gammaFluenceAbs[i]+=localRun->gammaFluenceAbs[i];
		gammaFluenceOver200Abs[i]+=localRun->gammaFluenceOver200Abs[i];
		gammaCreationAbs[i]+=localRun->gammaCreationAbs[i];
	}
	// mergeing pion location in the world
	const G4int nW = cfg_.nWorldX * cfg_.nWorldZ;
	for (int i=0;i<nW;i++)
	{
		pionFluenceWorld[i]+=localRun->pionFluenceWorld[i];
	}
	// mergeing spectra of pion in/out and gamma 
	for (int k=0;k<cfg_.energyBins;k++)
	{
		pionEnergyIn[k] += localRun->pionEnergyIn[k];
		gammaEnergy[k] += localRun->gammaEnergy[k];
		pionEnergyOut[k] += localRun->pionEnergyOut[k];
	
	}
	
	G4Run::Merge(aRun);
	
}
