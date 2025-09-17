#include "stepping.hh" 
#include "G4TrackStatus.hh"
namespace 
{
	constexpr G4int PDG_Gamma    = 22;     // γ
	constexpr G4int PDG_Electron = 11;     // e⁻
	constexpr G4int PDG_Positron = -11;    // e⁺
	constexpr G4int PDG_PionPlus = 211;    // π⁺
	constexpr G4int PDG_PionMinus= -211;   // π⁻
	constexpr G4int PDG_PionZero = 111;    // π⁰
	constexpr G4int PDG_Neutron  = 2112;   // n
	constexpr G4int maxEnergyBin = 1000; // 0-1 GeV, 1 MeV/bin
	constexpr int kMaxIdx = maxEnergyBin - 1; // 999
	constexpr G4double minAllowedStep = 1.1e-6 * mm;  // Prevents 0-length steps
	constexpr G4double minAllowedTime = 1.0 * microsecond; // Time limit for tracking particles
	auto clampIndex = [](int idx){return idx < 0 ? 0 : (idx > kMaxIdx ? kMaxIdx : idx);};
	static bool gRunPiPlusMain   = true;   // main physics: detect pi+ at end detector
	static bool gRunConvStats    = false;  // converter stats: gamma/e± flux, etc.
	static bool gRunDebug        = false;  // ad-hoc prints/checks
	inline bool IsAllowedPDG(G4int pdg)
    {
		// gamma, e-, e+, pi+, pi-, neutron
		switch (pdg) 
		{
			case PDG_Gamma: case PDG_Electron: case PDG_Positron: case PDG_PionPlus: case PDG_PionMinus: case PDG_Neutron: return true;
			default: return false;
		}
    }
	inline void gridIndex(const G4ThreeVector& p0, const G4ThreeVector& p1, int& i, int& j)
    {
		static const double invD = 1.0/(0.2*mm);
		const double xmid = 0.5*(p0.x()+p1.x());
		const double zmid = 0.5*(p0.z()+p1.z());
		i = int(xmid*invD + 50.0);
		j = int(zmid*invD + 100.0);
	}
}
MySteppingAction::~MySteppingAction()
{}

void MySteppingAction::UserSteppingAction(const G4Step *step)
{    
	auto* track = step->GetTrack();
	auto* def   = track->GetDefinition();
	const G4int    pdg    = def->GetPDGEncoding();
	
	if (!IsAllowedPDG(pdg)) { track->SetTrackStatus(fStopAndKill); return;	} // Kill any particle not in the allowed list
	
	const G4double ekin   = track->GetKineticEnergy();   
	const G4double gtime  = track->GetGlobalTime();
	
	if (step->GetStepLength() < minAllowedStep || gtime > minAllowedTime) { track->SetTrackStatus(fStopAndKill); return;} // Kill particles with too short steps or too long time

	if (pdg != PDG_PionPlus && ekin<200*MeV) { track->SetTrackStatus(fStopAndKill); return; } // Kill low-energy non-pion particles
    
	auto* pre = step->GetPreStepPoint();
	const G4VPhysicalVolume* prePV  = pre->GetPhysicalVolume();
	auto  IsFirstStepInVolume = step->IsFirstStepInVolume();
    
	if (pdg == PDG_PionPlus)
	{
		if (gRunPiPlusMain)
		{
			if (prePV==fAbsorberPV && IsFirstStepInVolume)
			{
				int i = clampIndex(static_cast<int>(ekin/MeV + 0.5));
				fRunAction->fRun->pionEnergyIn[i] +=1;
				fEventAction->ifPionPProduced=true;
							
			}
			G4ThreeVector startPos = pre->GetPosition();
			G4ThreeVector endPos = step->GetPostStepPoint()->GetPosition();
			int i,j;
			gridIndex(startPos, endPos, i, j);
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
				
			if (prePV==fWorldPV && IsFirstStepInVolume)
			{
				int i = clampIndex(static_cast<int>(ekin/MeV + 0.5));
				fRunAction->fRun->pionEnergyOut[i] +=1;
				fRunAction->fRun->npiPosOut+=1;
				
			}
		}
	}
    if (gRunConvStats)
	{
		auto  IsFirstStepInVolume = step->IsFirstStepInVolume();
		if (pdg == PDG_Gamma  && prePV==fAbsorberPV && IsFirstStepInVolume)
		{
	    	int i = clampIndex(static_cast<int>(ekin/MeV + 0.5));
	    	fRunAction->fRun->gammaEnergy[i] +=1;
		}

    	if (pdg == PDG_Electron&& prePV==fAbsorberPV)
    	{
			G4ThreeVector startPos = pre->GetPosition();
			G4ThreeVector endPos = step->GetPostStepPoint()->GetPosition();
			int i,j;
			gridIndex(startPos, endPos, i, j);
			if (j<200 && i<100 && i>=0 && j>=0)
			{
				fEventAction->efluxboolian[i][j] = 1;
			}
    	}
    
    
     	if (pdg == PDG_Gamma&& prePV==fAbsorberPV)
    	{
			G4ThreeVector startPos = pre->GetPosition();
			G4ThreeVector endPos = step->GetPostStepPoint()->GetPosition();
			int i,j;
			gridIndex(startPos, endPos, i, j);
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
    
	}
	if(gRunDebug)
	{	
		fRunAction->fRun->steps+=1; // Count all steps
	}
}
