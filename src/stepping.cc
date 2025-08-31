#include "stepping.hh" 
#include "G4TrackStatus.hh"

constexpr G4int PDG_Gamma    = 22;     // γ
constexpr G4int PDG_Electron = 11;     // e⁻
constexpr G4int PDG_Positron = -11;    // e⁺
constexpr G4int PDG_PionPlus = 211;    // π⁺
constexpr G4int PDG_PionMinus= -211;   // π⁻
constexpr G4int PDG_PionZero = 111;    // π⁰
constexpr G4int PDG_Neutron  = 2112;   // n

MySteppingAction::~MySteppingAction()
{}
namespace
{
    inline bool IsAllowedPDG(G4int pdg)
    {
		// gamma, e-, e+, pi+, pi-, neutron
		switch (pdg) 
		{
			case PDG_Gamma: case PDG_Electron: case PDG_Positron: case PDG_PionPlus: case PDG_PionMinus: case PDG_Neutron: return true;
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
	auto* pre = step->GetPreStepPoint();
	const G4VPhysicalVolume* prePV  = pre->GetPhysicalVolume();
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
    if (pdg != PDG_PionPlus) 
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
		&& prePV==fAbsorberPV
		&& IsFirstStepInVolume)
	    {
		fRunAction->fRun->nin+=1;
	    }
	    
	    // Count the number of gamma rays that undergo Compton scattering
	    // Requirements:
	    // (1) Gamma ray
	    // (2) Compton scattering
	    // (3) Inside the target material
	    if (pdg==PDG_Gamma && processName=="compt"
		&& prePV==fAbsorberPV)
	    {
		fRunAction->fRun->ncompton+=1;
	    }
	    if (pdg==PDG_Gamma && processName=="conv"
		&& prePV==fAbsorberPV)
	    {
		fRunAction->fRun->nconv+=1;
	    }
	    if (step->GetTrack()->GetTrackID() == 1 && pdg==PDG_Electron && IsFirstStepInVolume
		&& processName=="Transportation"
		&& prePV==fAbsorberPV)
	    {
		fRunAction->fRun->nTransportation+=1;
	    }
	    if(afterPhotonuclear)
	    {
	    	std::cout<< "particle: "               << pdg<< std::endl;
	    	afterPhotonuclear = false;
	    }
	    
	    if (step->GetTrack()->GetTrackID() == 1 && pdg==PDG_Gamma && IsFirstStepInVolume
		&& processName=="photonNuclear" 
		&& prePV==fAbsorberPV)
	    {
	    	//std::cout << "process: "                  << processName
	    	//    << std::endl;
	    	//afterPhotonuclear = true;
		fRunAction->fRun->nPhotoNuclear+=1;
	    }
	    if (step->GetTrack()->GetTrackID() == 1 && pdg==PDG_Gamma && IsFirstStepInVolume
		&& processName!="photonNuclear" && processName!="Transportation" && processName!="conv" && processName!="compt" 
		&& processName=="phot"					&& prePV==fAbsorberPV)
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
		
	    if (pdg == PDG_PionPlus || processName=="pi+Inelastic" )
	    {
	    	//std::cout<< "process: "                  << processName << std::endl;
	    	//std::cout << "pi+ energy out = " <<ekin <<" MeV" <<std::endl;
	    	fRunAction->fRun->npiPosInElas +=1;
	    
	    }
    	    if (pdg == PDG_PionPlus && prePV==fAbsorberPV && IsFirstStepInVolume)
	    {
	    	//std::cout << "got pion+" <<std::endl;
	    	//std::cout << "pi+ energy in = " <<ekin <<" MeV" <<std::endl;
	    	int i = ekin;
	    	fRunAction->fRun->pionEnergyIn[i] +=1;
	    	//fRunAction->fRun->npiPosIn+=1;
	    	fEventAction->ifPionPProduced=true;
	    	//std::cout << "kinetic energy" << step->GetTrack()->GetKineticEnergy() << std::endl;
	    	
	    }
	    if (pdg == PDG_PionZero && IsFirstStepInVolume
		&& prePV==fAbsorberPV)
	    {
	    	//std::cout << "got pion0" <<std::endl;
	    	//fRunAction->fRun->npiZerIn+=1;
	    	fEventAction->ifPionZProduced=true;
	    }
	    if (pdg == PDG_PionMinus && IsFirstStepInVolume
		&& prePV==fAbsorberPV)
	    {
	    	//std::cout << "got pion-" <<std::endl;
	    	
	    	
	    	int i = ekin;
	    	//std::cout << "pi- energy in = " <<ekin <<" MeV" <<std::endl << i+150 << std::endl;
	    	fRunAction->fRun->pionEnergyIn[i] +=1;
	    	//fRunAction->fRun->npiNegIn+=1;
	    	fEventAction->ifPionNProduced=true;
	    	
	    }
	    if (pdg == PDG_PionPlus ||pdg == PDG_PionMinus)
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
	    if (pdg == PDG_PionPlus  && prePV==fWorldPV && IsFirstStepInVolume)
	    {
	    	//std::cout << "pi+ energy out = " <<ekin <<" MeV" <<std::endl;
	    	int i = ekin;
	    	fRunAction->fRun->pionEnergyOut[i] +=1;
	    	//std::cout << "got pion+" <<std::endl;
	    	fRunAction->fRun->npiPosOut+=1;
	    	
	    }
	    if (pdg == PDG_PionZero && prePV==fWorldPV && IsFirstStepInVolume)
	    {
	    	//std::cout << "got pion0" <<std::endl;
	    	fRunAction->fRun->npiZerOut+=1;
	    	
	    }
	    if (pdg == PDG_PionMinus && prePV==fWorldPV && IsFirstStepInVolume)
	    {
	    	//std::cout << "pi- energy out = " <<ekin <<" MeV" <<std::endl;
	    	int i = ekin;
	    	fRunAction->fRun->pionEnergyOut[i] +=1;
	    	//std::cout << "got pion+" <<std::endl;
	    	fRunAction->fRun->npiNegOut+=1;
	    	
	    }
	    
    }
    
    /*if (pdg == PDG_Gamma  && prePV==fAbsorberPV && IsFirstStepInVolume)
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
    if (pdg == PDG_Electron&& sprePV==fAbsorberPV)
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
    
    
     if (pdg == PDG_Gamma&& prePV==fAbsorberPV)
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
