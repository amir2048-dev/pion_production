#include "run.hh"
#include "event.hh"
Run::Run()
{
    ncompton=0;
    nin=0;
    nTransportation=0;
    nconv=0;
    nPhotoNuclear=0;
    npiPosIn=0;
    npiNegIn=0;
    npiZerIn=0;
    npiPosOut=0;
    npiNegOut=0;
    npiZerOut=0;
    steps=0;
    nPhot=0;
    neventsOfPion=0;
    npiPosInElas =0;
    //G4int runID = run->GetRunID();
    //std::stringstream strRunID;
    //strRunID << runID;
    //G4String fileName = "Output" +strRunID.str() + ".csv";
}
//void Run::RecordEvent(const G4UserEventAction* event)
//{
//	const MyEventAction* myevent = static_cast<const MyEventAction*>(event);
//	G4int evtNb = myevent->GetEventID();
//	G4cout << myevent->ifPionPProduced << G4endl;
//	
//
//
//}
void Run::Merge(const G4Run* aRun)
{
	
	const Run* localRun = static_cast<const Run*>(aRun);
	// mergeing all counters
	nin +=localRun->nin;
	ncompton +=localRun->ncompton;
	nTransportation+=localRun->nTransportation;
	nconv+=localRun->nconv;
	nPhotoNuclear+=localRun->nPhotoNuclear;
	npiPosIn+=localRun->npiPosIn;
	npiNegIn+=localRun->npiNegIn;
	npiZerIn+=localRun->npiZerIn;
	npiPosOut+=localRun->npiPosOut;
	npiNegOut+=localRun->npiNegOut;
	npiZerOut+=localRun->npiZerOut;
	steps+=localRun->steps;
	nPhot+=localRun->nPhot;
	neventsOfPion+=localRun->neventsOfPion;
	npiPosInElas+=localRun->npiPosInElas;
	// mergeing pion/e/gamma location in the absorber
	for (int i=0;i<100;i++)
	{
		for (int j=0;j<200;j++)
		{
			pionflux[i][j]+=localRun->pionflux[i][j];
			eflux[i][j]+=localRun->eflux[i][j];
			gammaflux[i][j]+=localRun->gammaflux[i][j];
			gammacreation[i][j]+=localRun->gammacreation[i][j];
			gammafluxover200[i][j]+=localRun->gammafluxover200[i][j];
		}
	}
	// mergeing pion location in the entire area
	for (int i=0;i<500;i++)
	{
		for (int j=0;j<1000;j++)
		{
			pionfluxlarge[i][j]+=localRun->pionfluxlarge[i][j];
		}
	}
	// mergeing energy distribuition of pion in/out and gamma 
	for (int k=0;k<1000;k++)
	{
		pionEnergyIn[k] += localRun->pionEnergyIn[k];
		gammaEnergy[k] += localRun->gammaEnergy[k];
		pionEnergyOut[k] += localRun->pionEnergyOut[k];
	
	}
	
	G4Run::Merge(aRun);
	
}
