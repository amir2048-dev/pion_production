#include "runAction.hh"
#include <thread>
#include <fstream>

MyRunAction::MyRunAction()
{
    
}
MyRunAction::~MyRunAction()
{}

G4Run* MyRunAction::GenerateRun()
{
	fRun = new Run();
	return fRun;

}

void MyRunAction::BeginOfRunAction(const G4Run* run)
{
    
    //    man->OpenFile(fileName);
}
void MyRunAction::EndOfRunAction(const G4Run* run)
{
    const Run* simRun = static_cast<const Run*>(run);
    
    if (IsMaster())
    {
    
    

    G4cout<< "in: " << simRun->nin << G4endl <<"Transportation: "<<simRun->nTransportation<< G4endl<< "compton:" << simRun->ncompton << G4endl << "conv: " << simRun->nconv << G4endl << "Photoelectric: " << simRun->nPhot <<    			   G4endl << "PhotoNuclear: " << simRun->nPhotoNuclear << G4endl
    << "pi+ in absorber: "<< simRun->npiPosIn << G4endl << "pi- in absorber: "<< simRun->npiNegIn << G4endl << "pi0 in absorber: "<< simRun->npiZerIn << G4endl << "pion total in: "<< simRun->npiPosIn+simRun->npiNegIn+simRun->npiZerIn << G4endl 
    << "pi+ out absorber: "<< simRun->npiPosOut << G4endl << "pi- out of absorber: "<< simRun->npiNegOut << G4endl << "pi0 out of absorber: "<< simRun->npiZerOut << G4endl << "pion total out: "<< simRun->npiPosOut+simRun->npiNegOut+simRun->npiZerOut << G4endl
    << "events of pionproduction " << simRun->neventsOfPion << G4endl << "number of steps: " << simRun->steps << G4endl << "number of pi+ inelastic scatering" << simRun->npiPosInElas << G4endl;
    
    
    
    std::string folderName = "test"; 
    std::string csvFilePath_pionflux = folderName + "/pionflux.csv";
    int pionfluxold[100][200]={0};
    std::ifstream csvFile_pionflux(csvFilePath_pionflux);
    if (!csvFile_pionflux.is_open()) 
    {
	std::cerr << "Failed to open input.csv" << std::endl;
    }
    if (csvFile_pionflux.peek() != std::ifstream::traits_type::eof()) 
    {
	    
	    
	    std::string line;
	    int row = 0;
	    while (getline(csvFile_pionflux, line) && row < 100)
	     {
		std::stringstream ss(line);
		std::string cell;
		int col = 0;

		while (getline(ss, cell, ',') && col < 200)
		{
		    // Convert the string to an integer
		    pionfluxold[row][col] = std::stoi(cell);
		    col++;
		}
		
		row++;
	    }

	    // Close the file
	    csvFile_pionflux.close();
    }
    std::ofstream outputFile_pionflux(csvFilePath_pionflux);
    for(int i = 0; i < 100; ++i)
    	{
        	for(int j = 0; j < 200; ++j)
        	{
    			outputFile_pionflux << simRun->pionflux[i][j] + pionfluxold[i][j] << ",";
        	}
        	outputFile_pionflux << "\n";

        }   
    outputFile_pionflux.close();
    
    
    
 std::string csvFilePath_pionfluxlarge = folderName + "/pionfluxlarge.csv";
    int pionfluxlargeold[500][1000]={0};
    std::ifstream csvFile_pionfluxlarge(csvFilePath_pionfluxlarge);
    if (!csvFile_pionfluxlarge.is_open()) 
    {
	std::cerr << "Failed to open input.csv" << std::endl;
    }
    if (csvFile_pionfluxlarge.peek() != std::ifstream::traits_type::eof()) 
    {
	    
	    
	    std::string line;
	    int row = 0;
	    while (getline(csvFile_pionfluxlarge, line) && row < 500)
	     {
		std::stringstream ss(line);
		std::string cell;
		int col = 0;

		while (getline(ss, cell, ',') && col < 1000)
		{
		    // Convert the string to an integer
		    pionfluxlargeold[row][col] = std::stoi(cell);
		    col++;
		}
		
		row++;
	    }

	    // Close the file
	    csvFile_pionfluxlarge.close();
    }
    std::ofstream outputFile_pionfluxlarge(csvFilePath_pionfluxlarge);
    for(int i = 0; i < 500; ++i)
    	{
        	for(int j = 0; j < 1000; ++j)
        	{
    			outputFile_pionfluxlarge << simRun->pionfluxlarge[i][j] + pionfluxlargeold[i][j] << ",";
        	}
        	outputFile_pionfluxlarge << "\n";

        }   
    outputFile_pionfluxlarge.close();
     
    
    
    std::string csvFilePath_eflux = folderName + "/eflux.csv";
    int efluxold[100][200]={0};
    std::ifstream csvFile_eflux(csvFilePath_eflux);
    if (!csvFile_eflux.is_open()) 
    {
	std::cerr << "Failed to open input.csv" << std::endl;
    }
    if (csvFile_eflux.peek() != std::ifstream::traits_type::eof()) 
    {
	    
	    
	    std::string line;
	    int row = 0;
	    while (getline(csvFile_eflux, line) && row < 100)
	     {
		std::stringstream ss(line);
		std::string cell;
		int col = 0;

		while (getline(ss, cell, ',') && col < 200)
		{
		    // Convert the string to an integer
		    efluxold[row][col] = std::stoi(cell);
		    col++;
		}
		
		row++;
	    }

	    // Close the file
	    csvFile_eflux.close();
    }
    std::ofstream outputFile_eflux(csvFilePath_eflux);
    for(int i = 0; i < 100; ++i)
    	{
        	for(int j = 0; j < 200; ++j)
        	{
    			outputFile_eflux << simRun->eflux[i][j] + efluxold[i][j] << ",";
        	}
        	outputFile_eflux << "\n";

        }   
    outputFile_eflux.close();
    
    
    std::string csvFilePath_gammaflux = folderName + "/gammaflux.csv";
    int gammafluxold[100][200]={0};
    std::ifstream csvFile_gammaflux(csvFilePath_gammaflux);
    if (!csvFile_gammaflux.is_open()) 
    {
	std::cerr << "Failed to open input.csv" << std::endl;
    }
    if (csvFile_gammaflux.peek() != std::ifstream::traits_type::eof()) 
    {
	    
	    
	    std::string line;
	    int row = 0;
	    while (getline(csvFile_gammaflux, line) && row < 100)
	     {
		std::stringstream ss(line);
		std::string cell;
		int col = 0;

		while (getline(ss, cell, ',') && col < 200)
		{
		    // Convert the string to an integer
		    gammafluxold[row][col] = std::stoi(cell);
		    col++;
		}
		
		row++;
	    }

	    // Close the file
	    csvFile_gammaflux.close();
    }
    std::ofstream outputFile_gammaflux(csvFilePath_gammaflux);
    for(int i = 0; i < 100; ++i)
    	{
        	for(int j = 0; j < 200; ++j)
        	{
    			outputFile_gammaflux << simRun->gammaflux[i][j] + gammafluxold[i][j] << ",";
        	}
        	outputFile_gammaflux << "\n";

        }   
    outputFile_gammaflux.close();
    
    std::string csvFilePath_gammacreation = folderName + "/gammacreation.csv";
    int gammacreationold[100][200]={0};
    std::ifstream csvFile_gammacreation(csvFilePath_gammacreation);
    if (!csvFile_gammacreation.is_open()) 
    {
	std::cerr << "Failed to open input.csv" << std::endl;
    }
    if (csvFile_gammacreation.peek() != std::ifstream::traits_type::eof()) 
    {
	    
	    
	    std::string line;
	    int row = 0;
	    while (getline(csvFile_gammacreation, line) && row < 100)
	     {
		std::stringstream ss(line);
		std::string cell;
		int col = 0;

		while (getline(ss, cell, ',') && col < 200)
		{
		    // Convert the string to an integer
		    gammacreationold[row][col] = std::stoi(cell);
		    col++;
		}
		
		row++;
	    }

	    // Close the file
	    csvFile_gammacreation.close();
    }
    std::ofstream outputFile_gammacreation(csvFilePath_gammacreation);
    for(int i = 0; i < 100; ++i)
    	{
        	for(int j = 0; j < 200; ++j)
        	{
    			outputFile_gammacreation << simRun->gammacreation[i][j] + gammacreationold[i][j] << ",";
        	}
        	outputFile_gammacreation << "\n";

        }   
    outputFile_gammacreation.close();
    
    
    std::string csvFilePath_gammafluxover200 = folderName + "/gammafluxover200.csv";
    int gammafluxover200old[100][200]={0};
    std::ifstream csvFile_gammafluxover200(csvFilePath_gammafluxover200);
    if (!csvFile_gammafluxover200.is_open()) 
    {
	std::cerr << "Failed to open input.csv" << std::endl;
    }
    if (csvFile_gammafluxover200.peek() != std::ifstream::traits_type::eof()) 
    {
	    
	    
	    std::string line;
	    int row = 0;
	    while (getline(csvFile_gammafluxover200, line) && row < 100)
	     {
		std::stringstream ss(line);
		std::string cell;
		int col = 0;

		while (getline(ss, cell, ',') && col < 200)
		{
		    // Convert the string to an integer
		    gammafluxover200old[row][col] = std::stoi(cell);
		    col++;
		}
		
		row++;
	    }

	    // Close the file
	    csvFile_gammafluxover200.close();
    }
    std::ofstream outputFile_gammafluxover200(csvFilePath_gammafluxover200);
    for(int i = 0; i < 100; ++i)
    	{
        	for(int j = 0; j < 200; ++j)
        	{
    			outputFile_gammafluxover200 << simRun->gammafluxover200[i][j] + gammafluxover200old[i][j] << ",";
        	}
        	outputFile_gammafluxover200 << "\n";

        }   
    outputFile_gammafluxover200.close();
    
    
    std::string csvFilePath_in = folderName + "/pion_energy_in.csv";
    std::string csvFilePath_out = folderName + "/pion_energy_out.csv";
    std::string csvFilePath_gamma = folderName + "/gamma_energy_in.csv";
    int energy_in[1000]={0};
    int energy_out[1000]={0};
    int energy_gamma[1000]={0};
    G4cout<< "file opend" << G4endl;
    std::ifstream csvFile_in(csvFilePath_in);
    if (!csvFile_in.is_open()) 
    {
	std::cerr << "Failed to open input.csv" << std::endl;
    }
   
    if (csvFile_in.peek() != std::ifstream::traits_type::eof()) 
    {
    	std::string line;
	getline(csvFile_in, line);
	std::stringstream ss(line);
	std::string cell;
	int col = 0;

	while (getline(ss, cell, ',') && col < 1000)
	{
	    // Convert the string to an integer
	    energy_in[col] = std::stoi(cell);
	    col++;
	}
	    // Close the file
	    csvFile_in.close();
    }
    G4cout << "file cloesd" << G4endl;
    std::ifstream csvFile_out(csvFilePath_out);
     if (!csvFile_out.is_open()) 
    {
	std::cerr << "Failed to open input.csv" << std::endl;
    }
    if (csvFile_out.peek() != std::ifstream::traits_type::eof()) 
    {
    	std::string line;
	getline(csvFile_out, line);
	std::stringstream ss(line);
	std::string cell;
	int col = 0;

	while (getline(ss, cell, ',') && col < 1000)
	{
	    // Convert the string to an integer
	    energy_out[col] = std::stoi(cell);
	    col++;
	    
	}
	    // Close the file
	csvFile_out.close();
    }
    std::ifstream csvFile_gamma(csvFilePath_gamma);
     if (!csvFile_gamma.is_open()) 
    {
	std::cerr << "Failed to open gamma energy" << std::endl;
    }
    if (csvFile_gamma.peek() != std::ifstream::traits_type::eof()) 
    {
    	std::string line;
	getline(csvFile_gamma, line);
	std::stringstream ss(line);
	std::string cell;
	int col = 0;

	while (getline(ss, cell, ',') && col < 1000)
	{
	    // Convert the string to an integer
	    energy_gamma[col] = std::stoi(cell);
	    col++;
	    
	}
	    // Close the file
	csvFile_gamma.close();
    }
    std::ofstream outputFile_in(csvFilePath_in);
    for(int i = 0; i < 1000; ++i)
    	{
    		outputFile_in << simRun->pionEnergyIn[i] + energy_in[i] << ",";
    		
        }   
    outputFile_in.close();
    std::ofstream outputFile_out(csvFilePath_out);
    for(int i = 0; i < 1000; ++i)
    	{
    		outputFile_out << simRun->pionEnergyOut[i] + energy_out[i] << ",";
    		
        }   
    outputFile_out.close();
    std::ofstream outputFile_gamma(csvFilePath_gamma);
    for(int i = 0; i < 1000; ++i)
    	{
    		outputFile_gamma << simRun->gammaEnergy[i] + energy_gamma[i] << ",";
    		
        }   
    outputFile_out.close();
    
  
    G4cout
     << G4endl
     << "--------------------End of Global Run-----------------------"
     << G4endl
     << "  The run was " << simRun->nin << " events ";
  }
    //    // save histograms & ntuple
    //    //
    //    analysisManager->Write();
    //    analysisManager->CloseFile();
    //
    
}

