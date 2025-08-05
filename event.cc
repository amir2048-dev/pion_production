#include "event.hh"

MyEventAction::MyEventAction()
{
    fEdep = 0.;
    fEventID = 0;
    ifPionPProduced = false;
    ifPionNProduced = false;
    ifPionZProduced = false;
	
}
MyEventAction::~MyEventAction()
{
	//delete boolianArray;
}

void MyEventAction::BeginOfEventAction(const G4Event*)
{	
	if (ifPionPProduced || ifPionNProduced)
    	{
		for(int i = 0; i < 500; ++i)
    		{
       		 	for(int j = 0; j < 1000; ++j)
       		 	{
        			pionfluxlargeboolian[i][j] = 0; 
        		}
        	}  
        	for(int i = 0; i < 100; ++i)
    		{
			for(int j = 0; j < 200; ++j)
			{
		    		pionfluxboolian[i][j] = 0;
		    	}
        	} 
        	
        }
	fEdep = 0.;
	ifPionPProduced = false;
    	ifPionNProduced = false;
    	ifPionZProduced = false;
    
    	//IF I WANT TO TRACK ELECTRONS AND GAMMA USE THIS AND PUT THE FOR ABOVE IN COMENT
    	/*
    	for(int i = 0; i < 100; ++i)
    	{
        	for(int j = 0; j < 200; ++j)
        	{
            		pionfluxboolian[i][j] = 0;
            		efluxboolian[i][j] = 0;
            		gammafluxboolian[i][j] = 0;
            		gammafluxover200boolian[i][j] = 0;
            		gammacreationboolian[i][j] = 0;
            	}
        } 
        */
     
	
}

void MyEventAction::EndOfEventAction(const G4Event*)
{
	
	//G4cout<<fEdep<<G4endl;
//	G4AnalysisManager *man = G4AnalysisManager::Instance();
//	man->FillH1(0, fEdep);
//	man->FillNtupleDColumn(0, fEdep);
//	man->AddNtupleRow();
    Run* frun = static_cast <Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    fEventID++;
    if (ifPionPProduced || ifPionNProduced)
    {
    	if (ifPionPProduced)
    	{
    	frun->npiPosIn+=1;
    	}
    	if (ifPionNProduced)
    	{
    	frun->npiNegIn+=1;
    	}
    	frun->neventsOfPion++;
    	for(int i = 0; i < 100; ++i)
    	{
        	for(int j = 0; j < 200; ++j)
        	{
            		frun->pionflux[i][j] += pionfluxboolian[i][j];
            	}
        } 
        for(int i = 0; i < 500; ++i)
    	{
        	for(int j = 0; j < 1000; ++j)
        	{
            		frun->pionfluxlarge[i][j] += pionfluxlargeboolian[i][j];
            	}
        }     
    }
    if (ifPionZProduced)
    {
    	frun->npiZerIn+=1;
    }
    /*
    for(int i = 0; i < 100; ++i)
    	{
        	for(int j = 0; j < 200; ++j)
        	{
            		frun->eflux[i][j] += efluxboolian[i][j];
            		frun->gammaflux[i][j] += gammafluxboolian[i][j];
            		frun->gammafluxover200[i][j] += gammafluxover200boolian[i][j];
            		frun->gammacreation[i][j] += gammacreationboolian[i][j];
            	}
        }  
	*/
    //std::cout <<
    //std::endl << "\t\t End of event " << fEventID << ""
    //<< std::endl
    //<< "---------------------------------------------------------------------------" << //std::endl;
}
void MyEventAction::AddEdep(G4double edep)
{
	fEdep +=edep;
}

