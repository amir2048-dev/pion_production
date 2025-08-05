#include <iostream>
#include "G4RunManagerFactory.hh"
#include "G4MTRunManager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "construction.hh"
#include "physics.hh"
#include "PhysicsList.hh"
#include "action.hh"
#include "QBBC.hh"
#include "QGSP_BERT.hh"
#include "FTFP_BERT.hh"
#include "G4Timer.hh"
#include "Randomize.hh"
#include <ctime>

int main(int argc, char** argv)
{
    G4Timer myTimer;
    myTimer.Start(); 
    G4long seed = static_cast<G4long>(time(nullptr));
    CLHEP::HepRandom::setTheSeed(seed);
    
    auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
    //G4RunManagerFactory *runManager = new G4RunManagerFactory();
    runManager->SetUserInitialization(new MyDetectorConstruction());
    //runManager->SetUserInitialization(new MyphysicsList());
    
    G4VModularPhysicsList* physicsList = new PhysicsList;
    //G4VModularPhysicsList* physicsList = new QGSP_BERT;
    //G4VModularPhysicsList* physicsList = new FTFP_BERT;
    //	G4VModularPhysicsList* physicsList = new G4EmStandardPhysics;
    
    //physicsList->SetVerboseLevel(1);
    runManager->SetUserInitialization(physicsList);
    //G4StepLimiterPhysics* stepLimitPhys = new G4StepLimiterPhysics();
    //runManager->SetUserInitialization(stepLimitPhys);
    runManager->SetUserInitialization(new MyActionInitialization());
    G4UIExecutive *ui = 0;
    if (argc == 1)
    {
    	runManager->Initialize();
        ui = new G4UIExecutive(argc, argv);
    }
    G4VisExecutive *visManager = new G4VisExecutive();
    visManager->Initialize();
    G4UImanager *UImanager = G4UImanager::GetUIpointer();
    if (ui)
    {
        G4UImanager *UImanager = G4UImanager::GetUIpointer();
        UImanager->ApplyCommand("/control/execute vis.mac");
        ui->SessionStart();
    }
    else
    {
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command + fileName);
    }
    
    //UImanager->ApplyCommand("/control/cout/ignoreThreadsExcept 0");
    //UImanager->ApplyCommand("/event/verbose 0");
    //UImanager->ApplyCommand("/tracking/verbose 0");
    //UImanager->ApplyCommand("/run/verbose 2");
    //UImanager->ApplyCommand("/control/verbose 2");
    //UImanager->ApplyCommand("/run/printProgress 100");
    //UImanager->ApplyCommand("/run/beamOn 1000");
    
    //
    //job termination
    delete visManager;
    delete runManager;
    delete ui;
    
    myTimer.Stop();
    G4cout<<"Spend time:"<<myTimer<<G4endl;
    // return
    return 0;
    
}
