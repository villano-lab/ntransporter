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
#include <cmath>

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


typedef std::vector<G4double> doubles;

// thermal energy (25 C)
constexpr G4double Etherm = 0.025692579120652998*eV;


// kernel of a Maxwell-Boltzmann thermal distribution of energies
G4double MaxwellBoltzmannKernel(G4double E) {
  return std::sqrt(E)*std::exp(-E/Etherm);
}

// trapezoidal integration of (x,y) data
G4double trap(const doubles &x, const doubles &y) {
  size_t n = x.size();
  if (n != y.size()) {
    throw std::invalid_argument("Invalid arguments to function trap(). " 
      "Vectors must have the same length.");
  }
  G4double s = 0;
  for (size_t i = 1; i < n; ++i) {
    s += (y[i] + y[i-1])*(x[i] - x[i-1])/2.;
  }
  return s;
}


// initialize Geant data with configured geometry and detector
void initialize() {
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
}


int main(int argc, char *argv[]) {

  G4int G, ng;
  std::string output_file_base, material_name;

  // parse command line args
  if (argc < 3)  {
    throw std::runtime_error("Error in NT_XS: Not enough arguments."
         "\n    Usage: ./NT_XS output_file_base_path material " 
         "[ngroups=100] [points_per_group=10]");
  } else {
    try {
      output_file_base = argv[1];
      material_name = argv[2];
      G = 100; // number of groups
      ng = 10; // points per group
      if (argc > 3) {
        G = std::stoi(argv[3]);
        if (argc > 4) {
          ng = std::stoi(argv[4]);
        }
      }
    } catch (std::invalid_argument) {
      throw std::runtime_error("Error in NT_XS: Invalid arguments. ngroups " 
        "and points_per_group must be numeric."
        "\n    Usage: ./NT_XS output_file_base_path material "
        "[ngroups=100] [points_per_group=10]");
    }
  }

  
  std::ofstream G4cout_file("G4cout_redirected_output.txt");
  
  auto G4cout_oldbuf = G4cout.rdbuf();
  G4cout.rdbuf(G4cout_file.rdbuf());
  
  // bounds of fast energies
  const G4double Emin = 0.1*eV, Emax = 20.*MeV; 

  // alpha (common ratio of group boundaries)
  G4double alpha = std::pow(Emin/Emax, 1./G);


  // array of group boundaries (one thermal group)
  doubles Eg(G+2);
  Eg[0] = Emax;

  // calculate group boundaries
  for (G4int g = 1; g < G+1; ++g) {
    Eg[g] = Eg[g-1]*alpha;
  }
  // lower bound
  Eg[G+1] = 0.;
  

  std::cout << Eg[G] << " should equal " << Emin << std::endl;

  
  // vector of energies at which to evaluate approximated integrals
  doubles Eeval(ng+1);


  // vector of cross sections (first element always zero; xsec[g] refers to group g)
  doubles xsec(G+2);
  xsec[0] = 0.;


  //throw std::runtime_error("End of this bit");

  // configure physics processes
  initialize();

  // neutron singleton
  G4Neutron *theNeutron = G4Neutron::Definition();

  // pull material table from SuperSim
  const CDMSMaterialTable *theTable = CDMSMaterialTable::GetInstance();
  
  // pull material data
  G4Material *material = theTable->GetMaterial(material_name);
  if (!material) {
    throw std::invalid_argument("Error in NT_XS: invalid material. Could not find in CDMS or NIST material tables.");
  } else {
    std::cout << material_name << " should equal " << material->GetName() << std::endl;
    G4ElementVector *elmVector = material->GetElementVector();
    G4int nElm = (*elmVector).size();
    for (int i = 0; i < nElm; ++i) {
      std::cout << i+1 << ": " << (*elmVector)[i]->GetName() << std::endl;
    }
  }

  
  // get ProcessManager for the neutron
  G4ProcessManager *theMan = theNeutron->GetProcessManager();
  
  // get vector of neutron processes
  G4int nProc = theMan->GetProcessListLength();
  G4ProcessVector *processes = theMan->GetProcessList();

  // dynamic particle: set energy, momentum
  G4DynamicParticle *dynamicNeutron = new G4DynamicParticle(theNeutron, 
                                              G4ThreeVector(0.,0.,1.), 0.);

  // process info
  G4HadronicProcess *elasticProc = dynamic_cast<G4HadronicProcess*>(
      (*processes)[2]);
  G4HadronicProcess *inelasticProc = dynamic_cast<G4HadronicProcess*>(
      (*processes)[3]);
  G4HadronicProcess *captureProc = dynamic_cast<G4HadronicProcess*>(
      (*processes)[4]);

  if (!elasticProc || !inelasticProc || !captureProc) {
    throw std::runtime_error("Error: casting one or more processes as "
    "G4HadronicProcess failed. Exitting.");
  }

  std::cout << elasticProc->GetProcessName() << std::endl;
  std::cout << inelasticProc->GetProcessName() << std::endl;
  std::cout << captureProc->GetProcessName() << std::endl;
  

  //for (G4int i = 0; i < nProc; ++i) {
  //  // if process is a hadronic process, print cross section
  //  G4HadronicProcess *thisProc = dynamic_cast<G4HadronicProcess*>(
  //    (*processes)[i]);
  //  if (thisProc) {
  //    // calculate and print
  //    //x = thisProc->GetCrossSectionDataStore()
  //    //            ->GetCrossSection(dynamicNeutron, material);
  //    std::cout << "Process " << i << " is " << thisProc->GetProcessName() 
  //              << std::endl;
  //    //thisProc->ProcessDescription(std::cout);
  //  }
  //}

  std::cout << "Done" << std::endl;

  G4cout.rdbuf(G4cout_oldbuf);

  return 0;
}
