#include "stepping.hh" 
#include "G4TrackStatus.hh"

MySteppingAction::MySteppingAction(MyEventAction *eventAction, MyRunAction* runAction)
{
    fEventAction = eventAction;
    fRunAction = runAction;
    afterPhotonuclear=false;
    
    
}

MySteppingAction::~MySteppingAction()
{}
namespace
{
    inline bool IsAllowedPDG(G4int pdg)
    {
		// gamma, e-, e+, pi+, pi-, neutron
		switch (pdg) 
		{
			case 22: case 11: case -11: case 211: case -211: case 2112: return true;
			default: return false;
		}
    }
}
void MySteppingAction::UserSteppingAction(const G4Step *step)
{    
	auto* track = step->GetTrack();
	auto* def   = track->GetDefinition();
	const G4int    pdg    = def->GetPDGEncoding();
	const G4double ekin   = track->GetKineticEnergy();   
	const G4double gtime  = track->GetGlobalTime();
	if (!IsAllowedPDG(pdg)) 
	{
		// G4cout << "Kill PDG=" << pdg << G4endl;
		track->SetTrackStatus(fStopAndKill);
		return;
	}
	fRunAction->fRun->steps+=1;
	static const G4double minAllowedStep = 1.1e-6 * mm;  // Prevents 0-length steps
	if (step->GetStepLength() < minAllowedStep)
	{
		std::cout << "step length: " << step->GetStepLength() << std::endl;
	}
	if (step->GetStepLength() < minAllowedStep)
	{
		step->GetTrack()->SetTrackStatus(fStopAndKill);
		return;
	}
	G4double timeLimit = 1.0 * microsecond;
	if (track->GetGlobalTime() > timeLimit) 
	{
		track->SetTrackStatus(fStopAndKill);  // Kill the particle
		return;
	}
    if (pdg != 211) 
    {
    	
    	if (ekin<200*MeV)
    	{
    		track->SetTrackStatus(fStopAndKill);
    		return;
    	}
    	//GAMMA AND E+
    	else if (false)
    	{
    	    auto  IsFirstStepInVolume = step->IsFirstStepInVolume();
    	    const G4String processName = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
	    // Count the number of gamma rays that enter the target region
	    // Requirements:
	    // (1) Primary particle (gamma ray) = track ID == 1
	    // (2) Inside the target material
	    // (3) First time particle meets this volume
	    if (step->GetTrack()->GetTrackID() == 1
		&& step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "physAbsorber"
		&& IsFirstStepInVolume)
	    {
		fRunAction->fRun->nin+=1;
	    }
	    
	    // Count the number of gamma rays that undergo Compton scattering
	    // Requirements:
	    // (1) Gamma ray
	    // (2) Compton scattering
	    // (3) Inside the target material
	    if (pdg==22 && processName=="compt"
		&& step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "physAbsorber")
	    {
		fRunAction->fRun->ncompton+=1;
	    }
	    if (pdg==22 && processName=="conv"
		&& step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "physAbsorber")
	    {
		fRunAction->fRun->nconv+=1;
	    }
	    if (step->GetTrack()->GetTrackID() == 1 && pdg==11 && IsFirstStepInVolume
		&& processName=="Transportation"
		&& step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "physAbsorber")
	    {
		fRunAction->fRun->nTransportation+=1;
	    }
	    if(afterPhotonuclear)
	    {
	    	std::cout<< "particle: "               << pdg<< std::endl;
	    	afterPhotonuclear = false;
	    }
	    
	    if (step->GetTrack()->GetTrackID() == 1 && pdg==22 && IsFirstStepInVolume
		&& processName=="photonNuclear" 
		&& step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "physAbsorber")
	    {
	    	//std::cout << "process: "                  << processName
	    	//    << std::endl;
	    	//afterPhotonuclear = true;
		fRunAction->fRun->nPhotoNuclear+=1;
	    }
	    if (step->GetTrack()->GetTrackID() == 1 && pdg==22 && IsFirstStepInVolume
		&& processName!="photonNuclear" && processName!="Transportation" && processName!="conv" && processName!="compt" 
		&& processName=="phot"					&& step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "physAbsorber")
	    {
	    	fRunAction->fRun->nPhot+=1;
	    	//std::cout << "process: "                  << processName
	    	//    << std::endl;
	    }
	    
	    //G4cout << "ProcName: " << step->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName() << G4endl;
	    G4double edep = step->GetTotalEnergyDeposit();
	    fEventAction->AddEdep(edep);
    	}
    }
    else
    {
    	G4String name = track->GetParticleDefinition()->GetParticleName();
    	auto  IsFirstStepInVolume = step->IsFirstStepInVolume();
    	const G4String processName = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
	    if (pdg == 211 || processName=="pi+Inelastic" )
	    {
	    	//std::cout<< "process: "                  << processName << std::endl;
	    	//std::cout << "pi+ energy out = " <<ekin <<" MeV" <<std::endl;
	    	fRunAction->fRun->npiPosInElas +=1;
	    
	    }
    	    if (pdg == 211 && step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "physAbsorber"&& IsFirstStepInVolume)
	    {
	    	//std::cout << "got pion+" <<std::endl;
	    	//std::cout << "pi+ energy in = " <<ekin <<" MeV" <<std::endl;
	    	int i = ekin;
	    	fRunAction->fRun->pionEnergyIn[i] +=1;
	    	//fRunAction->fRun->npiPosIn+=1;
	    	fEventAction->ifPionPProduced=true;
	    	//std::cout << "kinetic energy" << step->GetTrack()->GetKineticEnergy() << std::endl;
	    	
	    }
	    if (pdg == 111 && IsFirstStepInVolume
		&& step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "physAbsorber")
	    {
	    	//std::cout << "got pion0" <<std::endl;
	    	//fRunAction->fRun->npiZerIn+=1;
	    	fEventAction->ifPionZProduced=true;
	    }
	    if (pdg == -221 && IsFirstStepInVolume
		&& step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "physAbsorber")
	    {
	    	//std::cout << "got pion-" <<std::endl;
	    	
	    	
	    	int i = ekin;
	    	//std::cout << "pi- energy in = " <<ekin <<" MeV" <<std::endl << i+150 << std::endl;
	    	fRunAction->fRun->pionEnergyIn[i] +=1;
	    	//fRunAction->fRun->npiNegIn+=1;
	    	fEventAction->ifPionNProduced=true;
	    	
	    }
	    if (pdg == 211 ||pdg == -221)
	    {
	    	//std::cout << "Step length: "              << step->GetStepLength()/CLHEP::centimeter << " cm"
		//    << std::endl;
		//std::cout << "Step pos: "  << step->GetPreStepPoint()->GetPosition()<< std::endl;
		G4ThreeVector startPos = step->GetPreStepPoint()->GetPosition();
		G4ThreeVector endPos = step->GetPostStepPoint()->GetPosition();
		double zstart = startPos[2];
		double zend = endPos[2];
		int j = (zstart+zend)/(0.2*mm)+100;
		double xstart = startPos[0];
		double xend = endPos[0];
		int i = (xstart+xend)/(0.2*mm)+50;
		//std::cout << "zstart " <<zstart  << "  zend " <<zend << "  index " << j << std::endl ;
		//std::cout << "xstart " <<xstart  << "  xend " <<xend << "  index " << i << std::endl ;
		if (j<200 && i<100 && i>=0 && j>=0)
		{
			fEventAction->pionfluxboolian[i][j] = 1;
		}
		if (false)
		{
			j +=399;
			i +=199; 
			fEventAction->pionfluxlargeboolian[i][j] = 1;
		}
	    	
	    }
	    if (pdg == 211  && step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "physWorld"&& IsFirstStepInVolume)
	    {
	    	//std::cout << "pi+ energy out = " <<ekin <<" MeV" <<std::endl;
	    	int i = ekin;
	    	fRunAction->fRun->pionEnergyOut[i] +=1;
	    	//std::cout << "got pion+" <<std::endl;
	    	fRunAction->fRun->npiPosOut+=1;
	    	
	    }
	    if (pdg == 111 && step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "physWorld"&& IsFirstStepInVolume)
	    {
	    	//std::cout << "got pion0" <<std::endl;
	    	fRunAction->fRun->npiZerOut+=1;
	    	
	    }
	    if (pdg == -221 && step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "physWorld"&& IsFirstStepInVolume)
	    {
	    	//std::cout << "pi- energy out = " <<ekin <<" MeV" <<std::endl;
	    	int i = ekin;
	    	fRunAction->fRun->pionEnergyOut[i] +=1;
	    	//std::cout << "got pion+" <<std::endl;
	    	fRunAction->fRun->npiNegOut+=1;
	    	
	    }
	    
    }
    
    /*if (pdg == 22  && step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "physAbsorber"&& IsFirstStepInVolume)
	{
	    	//std::cout << "e- energy out = " <<ekin <<" MeV" <<std::endl;
	    	int i = ekin;
	    	fRunAction->fRun->gammaEnergy[i] +=1;
	    	//std::cout << "got pion+" <<std::endl;
	//    	fRunAction->fRun->npiPosOut+=1;
	    	
	}*/
	/*
    if (false)
    {
    	    std::cout << "e- energy out = " <<ekin <<" MeV" <<std::endl;
	    std::cout
	    << "Track "                     << step->GetTrack()->GetTrackID()
	    << ", particle: "               << name
	    << std::endl
	    << "process: "                  << processName
	    << std::endl
	    << "pre step volume: "          << step->GetPreStepPoint()->GetPhysicalVolume()->GetName()
	    << std::endl
	    << "is first step in volume? "  << IsFirstStepInVolume
	    << std::endl
	    << "Step length: "              << step->GetStepLength()/CLHEP::centimeter << " cm"
	    << std::endl;
	}
   */
   /*
    if (pdg == 11&& step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "physAbsorber")
    {
    		
    		//std::cout << "e- Step length: "              << step->GetStepLength()/CLHEP::centimeter << " cm"
		//<< std::endl;
		//std::cout << "Step pos: "  << step->GetPreStepPoint()->GetPosition()<< std::endl;
		G4ThreeVector startPos = step->GetPreStepPoint()->GetPosition();
		G4ThreeVector endPos = step->GetPostStepPoint()->GetPosition();
		double zstart = startPos[2];
		double zend = endPos[2];
		int j = (zstart+zend)/(0.2*mm)+100;
		double xstart = startPos[0];
		double xend = endPos[0];
		int i = (xstart+xend)/(0.2*mm)+50;
		//std::cout << "zstart " <<zstart  << "  zend " <<zend << "  index " << j << std::endl ;
		//std::cout << "xstart " <<xstart  << "  xend " <<xend << "  index " << i << std::endl ;
		if (j<200 && i<100 && i>=0 && j>=0)
		{
			fEventAction->efluxboolian[i][j] = 1;
		}
    }
    
    
     if (pdg == 22&& step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "physAbsorber")
    {
    		//std::cout << "gamma Step length: "              << step->GetStepLength()/CLHEP::centimeter << " cm"
		//<< std::endl;
		//std::cout << "Step pos: "  << step->GetPreStepPoint()->GetPosition()<< std::endl;
		G4ThreeVector startPos = step->GetPreStepPoint()->GetPosition();
		G4ThreeVector endPos = step->GetPostStepPoint()->GetPosition();
		double zstart = startPos[2];
		double zend = endPos[2];
		int j = (zstart+zend)/(0.2*mm)+100;
		double xstart = startPos[0];
		double xend = endPos[0];
		int i = (xstart+xend)/(0.2*mm)+50;
		//std::cout << "zstart " <<zstart  << "  zend " <<zend << "  index " << j << std::endl ;
		//std::cout << "xstart " <<xstart  << "  xend " <<xend << "  index " << i << std::endl ;
		if (j<200 && i<100 && i>=0 && j>=0)
		{
			fEventAction->gammafluxboolian[i][j] = 1;
			if (ekin>200*MeV)
			{
				fEventAction->gammafluxover200boolian[i][j] = 1;
			}
			if (IsFirstStepInVolume)
			{
				fEventAction->gammacreationboolian[i][j] = 1;
			}
		}
    }
    */
    
    
    
}
