// initialize SuperSim with full physis list and print all process info for neutrons
//
// generate executable PROCINFO
//
// usage: ./PROCINFO

#include <iostream>

#include "SuperSim_Main.hh"
#include "G4UImanager.hh"

#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"

#include "G4Neutron.hh"

#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4VProcess.hh"
#include "G4HadronicProcess.hh"


#include "G4HadronElasticProcess.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4HadronFissionProcess.hh"

#include "G4ios.hh"

#include "NTUtilities.hh"


int main(int argc, char *argv[]) {

    std::ofstream G4cout_file("G4cout_procinfo_redirected_output.txt");

    auto G4cout_oldbuf = G4cout.rdbuf();
    G4cout.rdbuf(G4cout_file.rdbuf());

    // configure physics processes

    // manager object to configure everything
    std::cout << "Starting SuperSim_Main" << std::endl;
    SuperSim_Main *sMain = new SuperSim_Main();

    // set verbosity
    sMain->SetVerboseLevel(0);
    sMain->runManager.SetVerboseLevel(0);

    std::cout << "Configuring physics list" << std::endl;
    sMain->Configure("full"); // full Shielding physics list, etc

    // configure geometry and detector (required for initialization)
    std::cout << "Configuring environment geometry" << std::endl;
    sMain->UImanager->ApplyCommand("/CDMS/lab NoLab");
    sMain->UImanager->ApplyCommand("/CDMS/detector ZIP");

    // initialize: among other things, build physics processes and attach to 
    // particle process managers
    std::cout << "Initializing" << std::endl;
    sMain->runManager.Initialize();

    // neutron singleton
    std::cout << "Fetching neutron singleton" << std::endl;
    G4Neutron *theNeutron = G4Neutron::Definition();


    // get ProcessManager for the neutron
    std::cout << "Fetching neutron process manager" << std::endl;
    G4ProcessManager *theMan = theNeutron->GetProcessManager();

    // get vector of neutron processes
    G4ProcessVector *processes = theMan->GetProcessList();

    G4int nProc = theMan->GetProcessListLength();


    std::cout << nProc << " processes" << std::endl;
    for (G4int i = 0; i < nProc; ++i) {

        G4String procName = (*processes)[i]->GetProcessName();
        G4String procType = G4VProcess::GetProcessTypeName((*processes)[i]->GetProcessType());

        std::cout << "Process " << i << " : " << " name = " << procName << ", type = " << procType << std::endl;

        (*processes)[i]->DumpInfo();
        (*processes)[i]->ProcessDescription(std::cout);

    }


    delete sMain;
    G4cout.rdbuf(G4cout_oldbuf);

}


