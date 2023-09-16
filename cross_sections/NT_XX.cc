// calculate cross section group constants for isotropic simplified flux models for multiple materials/binnings
// 
// assume 1/E for fast neutrons, Maxwell-Boltzmann distribution at room temperature for thermal neutrons
// 
// generates executable NT_XX
// 
// usage:
// ./NT_XX output_file_base_path material1 material2 ... [-n ngroups1=100 ngroups2 ...]

#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <cmath>
#include <limits>

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

#include "NTUtilities.hh"


typedef std::vector<G4double> doubles;
typedef std::vector<G4int> ints;
typedef std::vector<std::string> strings;


// initialize Geant data with configured geometry and detector
void initialize() {
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
}


int main(int argc, char *argv[]) {

  G4int ng = 20;
  std::string output_file_base;
  strings material_names;
  ints Gs;

  // parse command line args
  if (argc < 3)  {
    throw std::runtime_error("Error in NT_XX: Not enough arguments."
         "\n    Usage: ./NT_XX output_file_base_path material1 material2 ... [-n ngroups1=100 ngroups2 ...]");
  } else {
    try {

      output_file_base = argv[1];

      size_t i = 2;
      while (i < argc && argv[i][0] != '-' && argv[i][1] != 'n') {
        material_names.push_back(argv[i]);
        ++i;
      }
      for (size_t j = i+1; j < argc; ++j) {
        Gs.push_back(std::stoi(argv[j]));
      }
      if (Gs.size() == 0) {
        Gs.push_back(100);
      }
    } catch (std::invalid_argument) {
      throw std::runtime_error("Error in NT_XX: Invalid arguments. ngroups must be numeric."
        "\n    Usage: ./NT_XX output_file_base_path material1 material2 ... [-n ngroups1=100 ngroups2 ...]");
    }
  }
  std::cout << std::endl << "material_names:" << std::endl;
  for (size_t i = 0; i < material_names.size(); ++i) {
    std::cout << material_names[i] << std::endl;
  }
  std::cout << std::endl << std::endl << "Ns" << std::endl;
  for (size_t i = 0; i < Gs.size(); ++i) {
    std::cout << Gs[i] << std::endl;
  }

  if (1) {
    return 0;
  }
  
  int G = Gs[0];
  std::string material_name = material_names[0];
  
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
  // lower bound (basically zero)
  Eg[G+1] = std::numeric_limits<G4double>::epsilon()*Eg[G];
  
  // vectors of energies, cross sections, and fluxes at which to evaluate approximated integrals
  doubles E_eval(ng+1);
  doubles xs_eval(ng+1); // scatters
  doubles xa_eval(ng+1); // absorption (captures)
  doubles phi_eval(ng+1);


  // vector of cross sections (first element always zero; xsec[g] refers to group g)
  doubles xs(G+2); // scatters
  xs[0] = 0.;
  doubles xt(G+2);
  xt[0] = 0.; // total


  // configure physics processes
  initialize();

  // neutron singleton
  std::cout << "Fetching neutron singleton" << std::endl;
  G4Neutron *theNeutron = G4Neutron::Definition();

  // pull material table from SuperSim
  std::cout << "Fetching material " << material_name << std::endl;
  const CDMSMaterialTable *theTable = CDMSMaterialTable::GetInstance();
  // pull material data
  G4Material *material = theTable->GetMaterial(material_name);
  if (!material) {
    throw std::invalid_argument("Error in NT_XX: invalid material \"" + material_name + "\". Could not find in CDMS or NIST material tables.");
  }

  
  // get ProcessManager for the neutron
  std::cout << "Fetching neutron process manager" << std::endl;
  G4ProcessManager *theMan = theNeutron->GetProcessManager();
  
  // get vector of neutron processes
  G4ProcessVector *processes = theMan->GetProcessList();

  // hadronic processes (except nFission, (*processes)[5])
  std::cout << "Fetching neutron hadronic processes" << std::endl;
  G4HadronicProcess *elasticProc = dynamic_cast<G4HadronicProcess*>(
      (*processes)[2]);
  G4HadronicProcess *inelasticProc = dynamic_cast<G4HadronicProcess*>(
      (*processes)[3]);
  G4HadronicProcess *captureProc = dynamic_cast<G4HadronicProcess*>(
      (*processes)[4]);

  if (!elasticProc || !inelasticProc || !captureProc) {
    throw std::runtime_error("Error in NT_XX: casting one or more processes as "
    "G4HadronicProcess failed. Exitting.");
  }

  G4CrossSectionDataStore *elasticDataStore = elasticProc->GetCrossSectionDataStore();
  G4CrossSectionDataStore *inelasticDataStore = inelasticProc->GetCrossSectionDataStore();
  G4CrossSectionDataStore *captureDataStore = captureProc->GetCrossSectionDataStore();

  // dynamic particle: set energy, momentum
  std::cout << "Creating dynamic neutron" << std::endl;
  G4DynamicParticle *dynamicNeutron = new G4DynamicParticle(theNeutron, 
                                              G4ThreeVector(0.,0.,1.), 0.);

  
  
  G4double group_min, group_max, r, phi_g;

  double (*phi_func)(double);
  
  std::cout << "Beginning group constant calculations" << std::endl;
  for (int g = 1; g < G+2; ++g) { // group g
    if (g == G+1) { // thermal group
      phi_func = MaxwellBoltzmannKernel;
    } else { // fast group
      phi_func = fast_kernel;
    }
    group_min = Eg[g]; // lower bound of group
    group_max = Eg[g-1]; // upper bound of group
    r = std::pow(group_max/group_min, 1./ng); // common ratio between evaluation points
    
    // evaluation points for integral
    E_eval[0] = group_min;
    for (int i = 1; i < ng; ++i) {
      E_eval[i] = E_eval[i-1]*r;
    }
    E_eval[ng] = group_max;

    // evaluate cross sections and fluxes
    for (int i = 0; i < E_eval.size(); ++i) {

      phi_eval[i] = phi_func(E_eval[i]);

      dynamicNeutron->SetKineticEnergy(E_eval[i]);

      xs_eval[i] = cm*phi_eval[i]
            *(elasticDataStore->GetCrossSection(dynamicNeutron, material) 
            + inelasticDataStore->GetCrossSection(dynamicNeutron, material));

      xa_eval[i] = cm*phi_eval[i]*captureDataStore->GetCrossSection(dynamicNeutron,
                                  material);

    }
    phi_g = trap(E_eval, phi_eval);
    
    xs[g] = trap(E_eval, xs_eval)/phi_g;
    xt[g] = xs[g] + trap(E_eval, xa_eval)/phi_g;
  }

  

  std::cout << "Group constants calculated" << std::endl;

  std::string filename = output_file_base + "_" 
                       + material_name + "_" 
                       + std::to_string(G) + "_" 
                       + std::to_string(ng) 
                       + "_xs.dat";

  std::cout << "Writing cross section data to " << filename << std::endl;

  std::ofstream outputStream(filename);

  outputStream << std::setprecision(17);

  for (int g = 0; g < G+2; ++g) {
    outputStream << g << " " << Eg[g] << " "<< xs[g] << " " << xt[g] << "\n";
  }

  outputStream.close();

  std::cout << "Done" << std::endl;

  G4cout.rdbuf(G4cout_oldbuf);

  return 0;
}
