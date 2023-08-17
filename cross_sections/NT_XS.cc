// calculate cross section group constants for isotropic simplified flux models
// 
// assume 1/E for fast neutrons, Maxwell-Boltzmann distribution at room 
// temperature for thermal neutrons
// 
// generates executable NT_XS
// 
// usage:
// ./NT_XS output_file_base_path material [ngroups=100] [points_per_group=10]

#include <iostream>
#include <vector>
#include <fstream>
#include <stdexcept>

#include "SuperSim_Main.hh"
#include "G4UImanager.hh"

#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"

#include "CDMSMaterialTable.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4ElementVector.hh"

#include "G4Neutron.hh"

#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4HadronicProcess.hh"

#include "G4ios.hh"


int main(int argc, char *argv[]) {

  G4int ngroups, points_per_group;
  std::string output_file;

  // parse command line args
  if (argc < 2)  {
    throw std::runtime_error("Error in NT_XS: Not enough arguments.\n    Usage: ./NT_XS output_file_base_path material [ngroups=100] [points_per_group=10]");
  } else {
    try {
      output_file = argv[1];
      ngroups = 100;
      points_per_group = 10;
      if (argc > 2) {
        ngroups = std::stoi(argv[2]);
        if (argc > 3) {
          points_per_group = std::stoi(argv[3]);
        }
      }
    } catch (std::invalid_argument) {
      throw std::runtime_error("Error in NT_XS: Invalid arguments. ngroups and points_per_group must be numeric.\n    Usage: ./NT_XS output_file_base_path material [ngroups=100] [points_per_group=10]");
    }
  }

  //std::cout << "output_file = " << output_file << std::endl << "ngroups = " << ngroups << std::endl << "points_per_group = " << points_per_group << std::endl;
  //throw std::runtime_error("End of this bit");
  
  std::ofstream file("G4cout_redirected_output.txt");
  
  auto G4cout_oldbuf = G4cout.rdbuf();
  G4cout.rdbuf(file.rdbuf());

  // manager object to configure everything
  SuperSim_Main *sMain = new SuperSim_Main();

  // set verbosity
  sMain->SetVerboseLevel(0);
  sMain->runManager.SetVerboseLevel(0);

  sMain->Configure("full"); // full Shielding physics list, etc
  
  // configure geometry and detector (required for initialization)
  sMain->UImanager->ApplyCommand("/CDMS/lab NoLab");
  sMain->UImanager->ApplyCommand("/CDMS/detector ZIP");
  
  // initialize: among other things, build physics processes and attach to 
  // particle process managers
  sMain->runManager.Initialize();

  // neutron singleton
  G4Neutron *theNeutron = G4Neutron::Definition();

  // pull material table from SuperSim
  const CDMSMaterialTable *theTable = CDMSMaterialTable::GetInstance();
  
  // pull material data
  G4Material *rock = theTable->GetMaterial("Norite");
  //G4int nElm = rock->GetNumberOfElements();
  //const G4ElementVector *elmVector = rock->GetElementVector();
  //const G4double *fracVector = rock->GetFractionVector();
  //const G4Element *elm = (*elmVector)[6]; // silicon

  
  // get ProcessManager for the neutron
  G4ProcessManager *theMan = theNeutron->GetProcessManager();
  
  // get vector of neutron processes
  G4int nProc = theMan->GetProcessListLength();
  G4ProcessVector *processes = theMan->GetProcessList();
  
  // variables to store cross sections, neutron energy
  G4double x, Eneut = 2.*keV;


  //std::cout << "Cross sections for " << elm->GetName() << " with " << Eneut 
  //          << " MeV neutron" << std::endl;
  
  // dynamic particle: set energy, momentum
  G4DynamicParticle *dynamicNeutron = new G4DynamicParticle(theNeutron, 
                                              G4ThreeVector(0.,0.,1.), Eneut);
  
  std::cout << nProc << " processes" << std::endl;
  for (G4int i = 0; i < nProc; ++i) {
    // if process is a hadronic process, print cross section
    G4HadronicProcess *thisProc = dynamic_cast<G4HadronicProcess*>((*processes)[i]);
    if (thisProc) {
      // calculate and print
      x = thisProc->GetCrossSectionDataStore()->GetCrossSection(dynamicNeutron, rock);
      std::cout << "Cross section for " << thisProc->GetProcessName() 
                << " : " << x << std::endl;
    }
  }

  std::cout << "Done" << std::endl;

  G4cout.rdbuf(G4cout_oldbuf);

  return 0;
}
