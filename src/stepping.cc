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
	inline void worldToIdx(double x, double z, double pixelX, double pixelZ, int nX, int nZ,  double xOrigin, double zOrigin, int& ix, int& iz)
	{
		const double gx = (x - xOrigin) / pixelX + 0.5 * nX;
		const double gz = (z - zOrigin) / pixelZ + 0.5 * nZ;
		ix = static_cast<int>(std::floor(gx));
		iz = static_cast<int>(std::floor(gz));
	}
	// Traverse all cells intersected by segment p0→p1 on X–Z plane.
	// Calls `visit(ix,iz)` for each crossed cell (no allocations).
	template <typename Visitor>
	inline int raycastCellsXZ(const G4ThreeVector& p0,const G4ThreeVector& p1, double pixelX, double pixelZ, int nX, int nZ, double xOrigin, double zOrigin, Visitor&& visit)
	{
		int ix0, iz0, ix1, iz1;
		worldToIdx(p0.x(), p0.z(), pixelX, pixelZ, nX, nZ, xOrigin, zOrigin, ix0, iz0);
		worldToIdx(p1.x(), p1.z(), pixelX, pixelZ, nX, nZ, xOrigin, zOrigin, ix1, iz1);

		if (ix0 == ix1 && iz0 == iz1) { visit(ix0, iz0); return 1; }

		const double dx = p1.x() - p0.x();
		const double dz = p1.z() - p0.z();
		const int stepx = (dx > 0) - (dx < 0);
		const int stepz = (dz > 0) - (dz < 0);

		const double invDx = (dx != 0.0) ? (1.0 / dx) : 1e300;
		const double invDz = (dz != 0.0) ? (1.0 / dz) : 1e300;

		const double xEdge0 = xOrigin + ( (ix0 + (stepx>0)) - 0.5*nX ) * pixelX;
		const double zEdge0 = zOrigin + ( (iz0 + (stepz>0)) - 0.5*nZ ) * pixelZ;

		double tMaxX   = (dx != 0.0) ? ( (xEdge0 - p0.x()) * invDx ) : 1e300;
		double tMaxZ   = (dz != 0.0) ? ( (zEdge0 - p0.z()) * invDz ) : 1e300;
		const double tDeltaX = (dx != 0.0) ? ( pixelX * std::abs(invDx) ) : 1e300;
		const double tDeltaZ = (dz != 0.0) ? ( pixelZ * std::abs(invDz) ) : 1e300;

		int ix = ix0, iz = iz0, count = 0;
		visit(ix, iz); ++count;

		// March until we reach end cell
		while (ix != ix1 || iz != iz1) {
			if (tMaxX < tMaxZ) { ix += stepx; tMaxX += tDeltaX; }
			else               { iz += stepz; tMaxZ += tDeltaZ; }
			visit(ix, iz); ++count;

			if (count > (nX + nZ + 8)) break; // safety
		}
		return count;
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
			G4ThreeVector p0 = pre->GetPosition();
			G4ThreeVector p1 = step->GetPostStepPoint()->GetPosition();
			raycastCellsXZ(p0, p1, cfg_.pixelX, cfg_.pixelZ, cfg_.nAbsorberX, cfg_.nAbsorberZ,
               cfg_.absorberXOrigin, cfg_.absorberZOrigin,
               [&](int ix, int iz){
                 if (0<=ix && ix<cfg_.nAbsorberX && 0<=iz && iz<cfg_.nAbsorberZ)
                   fEventAction->pionfluxboolian[ix][iz] = 1;
               });
			if (false)
			{
				raycastCellsXZ(p0, p1, cfg_.pixelX, cfg_.pixelZ, cfg_.nWorldX, cfg_.nWorldZ,
               cfg_.worldXOrigin, cfg_.worldZOrigin,
               [&](int ix, int iz){
                 if (0<=ix && ix<cfg_.nWorldX && 0<=iz && iz<cfg_.nWorldZ)
                   fEventAction->pionfluxlargeboolian[ix][iz] = 1;
               });
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
			G4ThreeVector p0 = pre->GetPosition();
			G4ThreeVector p1 = step->GetPostStepPoint()->GetPosition();
			raycastCellsXZ(p0, p1, cfg_.pixelX, cfg_.pixelZ, cfg_.nAbsorberX, cfg_.nAbsorberZ,
               cfg_.absorberXOrigin, cfg_.absorberZOrigin,
               [&](int ix, int iz){
                 if (0<=ix && ix<cfg_.nAbsorberX && 0<=iz && iz<cfg_.nAbsorberZ)
                   fEventAction->efluxboolian[ix][iz] = 1;
               });
    	}
    
    
     	if (pdg == PDG_Gamma&& prePV==fAbsorberPV)
    	{
			G4ThreeVector p0 = pre->GetPosition();
			G4ThreeVector p1 = step->GetPostStepPoint()->GetPosition();
			if (ekin>200*MeV)
			{
				raycastCellsXZ(p0, p1, cfg_.pixelX, cfg_.pixelZ, cfg_.nAbsorberX, cfg_.nAbsorberZ, cfg_.absorberXOrigin, cfg_.absorberZOrigin,
				[&](int ix, int iz){
				if (0<=ix && ix<cfg_.nAbsorberX && 0<=iz && iz<cfg_.nAbsorberZ)
				fEventAction->gammafluxover200boolian[ix][iz] = 1;
			});
			}
			if (IsFirstStepInVolume)
			{
				raycastCellsXZ(p0, p1, cfg_.pixelX, cfg_.pixelZ, cfg_.nAbsorberX, cfg_.nAbsorberZ, cfg_.absorberXOrigin, cfg_.absorberZOrigin,
				[&](int ix, int iz){
				if (0<=ix && ix<cfg_.nAbsorberX && 0<=iz && iz<cfg_.nAbsorberZ)
				fEventAction->gammacreationboolian[ix][iz] = 1;
				});
			}
			
    	}
    
	}

	if(gRunDebug)
	{	
		fRunAction->fRun->steps+=1; // Count all steps
	}
}
