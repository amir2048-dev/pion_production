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
#include "SimConfig.hh"

int main(int argc, char** argv)
{
    G4Timer myTimer;
    myTimer.Start(); 
    G4long seed = static_cast<G4long>(time(nullptr));
    CLHEP::HepRandom::setTheSeed(seed);
    SimConfig cfg;
    std::string macroPath;
    auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
    runManager->SetUserInitialization(new MyDetectorConstruction(cfg));
    G4VModularPhysicsList* physicsList = new PhysicsList;
    runManager->SetUserInitialization(physicsList);
    
    G4UIExecutive *ui = 0;
    if (argc == 1)
    {
        ui = new G4UIExecutive(argc, argv);
        macroPath = "vis.mac";
    }
    else
    {
        macroPath = argv[1];
    }
    runManager->SetUserInitialization(new MyActionInitialization(cfg,macroPath));
    
    
    G4VisExecutive *visManager = new G4VisExecutive();
    visManager->Initialize();
    G4UImanager *UImanager = G4UImanager::GetUIpointer();
    if (ui)
    {
        runManager->Initialize();
        UImanager->ApplyCommand("/control/execute vis.mac");
        ui->SessionStart();
    }
    else
    {
        G4String command = "/control/execute ";
        UImanager->ApplyCommand(command + macroPath);
    }
    
    
    //job termination
    delete visManager;
    delete runManager;
    delete ui;
    
    myTimer.Stop();
    G4cout<<"Spend time:"<<myTimer<<G4endl;
    // return
    return 0;
    
}
