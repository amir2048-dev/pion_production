#include "stepping.hh" 
#include "G4TrackStatus.hh"
#include <cmath>

namespace 
{
	auto clampIndex = [](int idx, const SimConfig& cfg){return idx < 0 ? 0 : (idx > cfg.energymaxIndex ? cfg.energymaxIndex : idx);};
	auto energyToIndex = [](G4double ekin, const SimConfig& cfg) -> int {
		return static_cast<int>(ekin / cfg.getEnergyBinWidth() + 0.5);
	};
	// Traverse segment p0->p1 across an X–Z grid. For each crossed cell (ix,iz)
	// calls visit(ix, iz, segLenInCell) where segLenInCell is the physical length
	// of the segment INSIDE that cell (same units as p0/p1, typically mm).
	template <typename Visitor>
	inline int raycastAccumulateXZ(const G4ThreeVector& p0,
								const G4ThreeVector& p1,
								double pixelX, double pixelZ,
								int nX, int nZ,
								double xOrigin, double zOrigin,
								Visitor&& visit)
	{
		// Full segment
		const double x0 = p0.x(), z0 = p0.z();
		const double x1 = p1.x(), z1 = p1.z();
		const double dx = x1 - x0, dz = z1 - z0;

		// Handle degenerate segment
		const double segLen = std::sqrt(dx*dx + dz*dz);
		if (segLen == 0.0) return 0;

		// Helper to map world → index
		auto w2i = [&](double x, double z, int& ix, int& iz)
		{
			const double gx = (x - xOrigin) / pixelX + 0.5 * nX;
			const double gz = (z - zOrigin) / pixelZ + 0.5 * nZ;
			ix = static_cast<int>(std::floor(gx));
			iz = static_cast<int>(std::floor(gz));
		};

		int ix, iz, ix1, iz1;
		w2i(x0, z0, ix, iz);
		w2i(x1, z1, ix1, iz1);

		// Fast path: segment stays within one cell
		if (ix == ix1 && iz == iz1)
		{
			if (0<=ix && ix<nX && 0<=iz && iz<nZ) visit(ix, iz, segLen);
			return 1;
		}

		// DDA setup
		const int stepx = (dx > 0) - (dx < 0);
		const int stepz = (dz > 0) - (dz < 0);
		const double invDx = (dx != 0.0) ? (1.0 / dx) : 1e300;
		const double invDz = (dz != 0.0) ? (1.0 / dz) : 1e300;

		// Position of the next grid boundary in x and z directions
		const double xEdge0 = xOrigin + ((ix + (stepx>0)) - 0.5* nX) * pixelX;
		const double zEdge0 = zOrigin + ((iz + (stepz>0)) - 0.5* nZ) * pixelZ;

		double tMaxX   = (dx != 0.0) ? ((xEdge0 - x0) * invDx) : 1e300;
		double tMaxZ   = (dz != 0.0) ? ((zEdge0 - z0) * invDz) : 1e300;
		const double tDeltaX = (dx != 0.0) ? ( pixelX * std::abs(invDx) ) : 1e300;
		const double tDeltaZ = (dz != 0.0) ? ( pixelZ * std::abs(invDz) ) : 1e300;

		// March along the segment; t is the param along p(t) = p0 + t*(p1-p0)
		double tPrev = 0.0;
		int count = 0;

		// We'll iterate until we pass t=1.0
		while (tPrev < 1.0)
		{
			// Decide which grid boundary we hit next
			double tNext;
			if (tMaxX < tMaxZ)
			{
				tNext = std::min(tMaxX, 1.0);
				const double lenCell = (tNext - tPrev) * segLen;
				if (0<=ix && ix<nX && 0<=iz && iz<nZ && lenCell > 0.0) visit(ix, iz, lenCell);
				tPrev = tNext;
				ix += stepx;
				tMaxX += tDeltaX;
			} 
			else 
			{
				tNext = std::min(tMaxZ, 1.0);
				const double lenCell = (tNext - tPrev) * segLen;
				if (0<=ix && ix<nX && 0<=iz && iz<nZ && lenCell > 0.0) visit(ix, iz, lenCell);
				tPrev = tNext;
				iz += stepz;
				tMaxZ += tDeltaZ;
			}
			++count;
			if (count > (nX + nZ + 8)) break; // safety against degenerate cases
		}

		return count;
	}
	// helpers for angle histogram binning using small-angle approximation
	// theta_x and theta_y are the angles in radians from the z-axis projected onto x-z and y-z planes
	inline int thetaXBin(double thetaX, double minThetaX, double maxThetaX, int nBins)
	{
		if (thetaX < minThetaX || thetaX >= maxThetaX) return -1;
		const double dThetaX = (maxThetaX - minThetaX) / nBins;
		int idx = static_cast<int>((thetaX - minThetaX) / dThetaX);
		if (idx < 0) idx = 0;
		else if (idx >= nBins) idx = nBins - 1;
		return idx;
	}
	inline int thetaYBin(double thetaY, double minThetaY, double maxThetaY, int nBins)
	{
		if (thetaY < minThetaY || thetaY >= maxThetaY) return -1;
		const double dThetaY = (maxThetaY - minThetaY) / nBins;
		int idx = static_cast<int>((thetaY - minThetaY) / dThetaY);
		if (idx < 0) idx = 0;
		else if (idx >= nBins) idx = nBins - 1;
		return idx;
	}
}

MySteppingAction::~MySteppingAction()
{}

void MySteppingAction::UserSteppingAction(const G4Step *step)
{    
	auto* track = step->GetTrack();
	auto* def   = track->GetDefinition();
	const G4int    pdg    = def->GetPDGEncoding();
	
	// Kill particles not in the allowed PDG list (if enabled)
	if (cfg_.killDisallowedPDG && !IS_ALLOWED_PDG(pdg)) { 
		track->SetTrackStatus(fStopAndKill); 
		return;
	}
	
	const G4double ekin   = track->GetKineticEnergy();   
	const G4double gtime  = track->GetGlobalTime();
	
	// Kill particles exceeding maximum time (if enabled)
	if (cfg_.killByTime && gtime > cfg_.maxTime) { 
		track->SetTrackStatus(fStopAndKill); 
		return;
	}

	// Kill low-energy non-pion particles (if enabled)
	if (cfg_.killLowEnergyNonPion && pdg != AllowedParticles::PionPlus && ekin < cfg_.lowEnergyThreshold) { 
		track->SetTrackStatus(fStopAndKill); 
		return;
	}
    
	auto* pre = step->GetPreStepPoint();
	const G4VPhysicalVolume* prePV  = pre->GetPhysicalVolume();
	auto* post = step->GetPostStepPoint();
	const G4VPhysicalVolume* postPV = post->GetPhysicalVolume();
	auto  IsFirstStepInVolume = step->IsFirstStepInVolume();
    
	if (pdg == AllowedParticles::PionPlus)
	{
		if (cfg_.runPiPlusMain)
		{
			if (prePV==fAbsorberPV && IsFirstStepInVolume)
			{
				int i = clampIndex(energyToIndex(ekin, cfg_), cfg_);
				fRunAction->fRun->pionEnergyIn[i] +=1;
				fRunAction->fRun->npiPosIn+=1;
							
			}
			G4ThreeVector p0 = pre->GetPosition();
			G4ThreeVector p1 = step->GetPostStepPoint()->GetPosition();
			raycastAccumulateXZ(p0, p1, cfg_.pixelX, cfg_.pixelZ, cfg_.nAbsorberX, cfg_.nAbsorberZ,
               cfg_.absorberXOrigin, cfg_.absorberZOrigin,
               [&](int ix, int iz, double lenCell){
                 if (0<=ix && ix<cfg_.nAbsorberX && 0<=iz && iz<cfg_.nAbsorberZ)
				 {
				   const int k = fRunAction->fRun->AbsIndex(ix, iz);
                   fRunAction->fRun->pionFluenceAbs[k] += lenCell;
				 }
               });
			if (cfg_.runWorldMap)
			{
				raycastAccumulateXZ(p0, p1, cfg_.worldPixelX, cfg_.worldPixelZ, cfg_.nWorldX, cfg_.nWorldZ,
               cfg_.worldXOrigin, cfg_.worldZOrigin,
               [&](int ix, int iz, double lenCell){
                 if (0<=ix && ix<cfg_.nWorldX && 0<=iz && iz<cfg_.nWorldZ)
				 {
                   const int k = fRunAction->fRun->WorldIndex(ix, iz);
                   fRunAction->fRun->pionFluenceWorld[k] += lenCell;
				 }
               });
			}
			if (prePV==fAbsorberPV && postPV!=fAbsorberPV)
			{
				int i = clampIndex(energyToIndex(ekin, cfg_), cfg_);
				fRunAction->fRun->pionEnergyOut[i] +=1;
				fRunAction->fRun->npiPosOut+=1;
			}
		}
	}
	
	// Collect theta = arccos(p.z / |p|) for particles leaving absorber (forward hemisphere only)
	if (prePV==fAbsorberPV && postPV!=fAbsorberPV && IS_ALLOWED_PDG(pdg))
	{
		G4ThreeVector momentum = track->GetMomentum();
		if (momentum.z() > 0.0) {
			double p_mag = momentum.mag();
			if (p_mag > 0.0) {
				double theta = std::acos(std::min(1.0, momentum.z() / p_mag)); // clamp to avoid numerical errors
				fRunAction->fRun->particleThetaHistograms[pdg].push_back(theta);
			}
		}
	}
	
	// Exit plane - store position and momentum for post-processing
	if (cfg_.enableExitPlane && prePV==fExitPlanePV)
	{
		if (pdg == AllowedParticles::PionPlus || cfg_.angleIncludeBackground)
		{
			G4ThreeVector momentum = track->GetMomentum();
			// Only consider forward-going tracks (pz > 0)
			if (momentum.z() > 0.0) {
				G4ThreeVector position = pre->GetPosition();
				G4ThreeVector p_norm = momentum.unit();
				Run::ExitPlaneHit hit{position.x(), position.y(), p_norm.x(), p_norm.y(), p_norm.z()};
				
				if (pdg == AllowedParticles::PionPlus)
				{
					fRunAction->fRun->pionExitPlanePositions.push_back(hit);
				}
				else
				{
					fRunAction->fRun->backgroundExitPlanePositions.push_back(hit);
				}
			}
		}
	}
    if (cfg_.runConvStats)
	{
		auto  IsFirstStepInVolume = step->IsFirstStepInVolume();
		if (pdg == AllowedParticles::Gamma  && prePV==fAbsorberPV && IsFirstStepInVolume)
		{
	    	int i = clampIndex(energyToIndex(ekin, cfg_), cfg_);
	    	fRunAction->fRun->gammaEnergy[i] +=1;
		}

    	if (pdg == AllowedParticles::Electron&& prePV==fAbsorberPV)
    	{
			G4ThreeVector p0 = pre->GetPosition();
			G4ThreeVector p1 = step->GetPostStepPoint()->GetPosition();
			raycastAccumulateXZ(p0, p1, cfg_.pixelX, cfg_.pixelZ, cfg_.nAbsorberX, cfg_.nAbsorberZ,
               cfg_.absorberXOrigin, cfg_.absorberZOrigin,
               [&](int ix, int iz, double lenCell){
                 if (0<=ix && ix<cfg_.nAbsorberX && 0<=iz && iz<cfg_.nAbsorberZ)
				 {
                   const int k = fRunAction->fRun->AbsIndex(ix, iz);
                   fRunAction->fRun->eFluenceAbs[k] += lenCell;
				 }
               });
    	}
    
    
     	if (pdg == AllowedParticles::Gamma&& prePV==fAbsorberPV)
    	{
			G4ThreeVector p0 = pre->GetPosition();
			G4ThreeVector p1 = step->GetPostStepPoint()->GetPosition();
			raycastAccumulateXZ(p0, p1, cfg_.pixelX, cfg_.pixelZ, cfg_.nAbsorberX, cfg_.nAbsorberZ, cfg_.absorberXOrigin, cfg_.absorberZOrigin,
				[&](int ix, int iz, double lenCell){
				if (0<=ix && ix<cfg_.nAbsorberX && 0<=iz && iz<cfg_.nAbsorberZ)
				{
					const int k = fRunAction->fRun->AbsIndex(ix, iz);
					fRunAction->fRun->gammaFluenceAbs[k] += lenCell;
				}
				});
			if (ekin>200*MeV)
			{
				raycastAccumulateXZ(p0, p1, cfg_.pixelX, cfg_.pixelZ, cfg_.nAbsorberX, cfg_.nAbsorberZ, cfg_.absorberXOrigin, cfg_.absorberZOrigin,
				[&](int ix, int iz, double lenCell){
				if (0<=ix && ix<cfg_.nAbsorberX && 0<=iz && iz<cfg_.nAbsorberZ)
				{
					const int k = fRunAction->fRun->AbsIndex(ix, iz);
					fRunAction->fRun->gammaFluenceOver200Abs[k] += lenCell;
				}
				});
			}
			if (IsFirstStepInVolume)
			{
				raycastAccumulateXZ(p0, p1, cfg_.pixelX, cfg_.pixelZ, cfg_.nAbsorberX, cfg_.nAbsorberZ, cfg_.absorberXOrigin, cfg_.absorberZOrigin,
				[&](int ix, int iz, double lenCell){
				if (0<=ix && ix<cfg_.nAbsorberX && 0<=iz && iz<cfg_.nAbsorberZ)
				{
					const int k = fRunAction->fRun->AbsIndex(ix, iz);
					fRunAction->fRun->gammaCreationAbs[k] += lenCell;
				}
				});
			}
			
    	}
    
	}
	
	if(cfg_.runDebug)
	{	
		fRunAction->fRun->steps+=1; // Count all steps
		fRunAction->fRun->debugFeature+=1; // Count all steps for debug feature
	}
}
