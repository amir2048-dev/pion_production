#include "generator.hh"
#include "Randomize.hh"

MyPrimaryGenerator::MyPrimaryGenerator(const SimConfig& cfg): cfg_(cfg)
{
	fParticleGun = new G4ParticleGun;
	G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName = "e-";
	G4ParticleDefinition *particle = particleTable->FindParticle(particleName);
	G4ThreeVector mom(0.,0.,1.);
	fParticleGun->SetParticleMomentumDirection(mom);
	fParticleGun->SetParticleDefinition(particle);
}

MyPrimaryGenerator::~MyPrimaryGenerator()
{
	delete fParticleGun; 
}

void MyPrimaryGenerator::GeneratePrimaries(G4Event *anEvent)
{
	G4double fact=1;
	G4double boxhalfsize = 1*cm;
	G4double beamradius = 0.3*cm;
	G4double z0 = -boxhalfsize; //in absorber
	G4double x0 = 0;
	G4double y0 = 0;
	
	while(true)
	{
		x0 = (G4UniformRand()-0.5)*beamradius*2;
  		y0 = (G4UniformRand()-0.5)*beamradius*2;
  		if ((x0*x0+y0*y0)<=beamradius*beamradius)
  		{
  			break;
  		}	
  	}
  	
  	G4double E = 0;
  	G4double y = 0;
  	G4double fE = 0;
  	G4double t = 147.28;
  	G4double ymax = std::pow(t/2,0.5)*std::pow(t,-1.5)*std::exp(-1/2);
  	
  
  	while (true)
  	{
  		E = G4UniformRand()*GeV;
  		y = G4UniformRand();
  		fE = (std::pow(E,0.5)*std::pow(t,-1.5)*std::exp(-E/t))/ymax;
  		if (y<fE && E>0.2*GeV)
  		{
 			break; 		
  		}
   	}
   	
   	//E=300*MeV; //CHANGE!
  	fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
	fParticleGun->GeneratePrimaryVertex(anEvent);
	fParticleGun->SetParticleEnergy(E);
}
