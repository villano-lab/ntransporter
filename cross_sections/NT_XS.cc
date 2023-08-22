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
    throw std::invalid_argument("Error in NT_XS: invalid material \"" + material_name + "\". Could not find in CDMS or NIST material tables.");
  }

  
  // get ProcessManager for the neutron
  std::cout << "Fetching neutron process manager" << std::cout;
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
    throw std::runtime_error("Error in NT_XS: casting one or more processes as "
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
  
  std::cout << "Beginning fast group constant calculations" << std::endl;
  for (int g = 1; g < G+1; ++g) { // group g (fast groups)
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

      phi_eval[i] = 1/E_eval[i];

      dynamicNeutron->SetKineticEnergy(E_eval[i]);

      xs_eval[i] = phi_eval[i]
            *(elasticDataStore->GetCrossSection(dynamicNeutron, material) 
            + inelasticDataStore->GetCrossSection(dynamicNeutron, material));

      xa_eval[i] = phi_eval[i]*captureDataStore->GetCrossSection(dynamicNeutron,
                                  material);

    }
    phi_g = trap(E_eval, phi_eval);
    
    xs[g] = trap(E_eval, xs_eval)/phi_g;
    xt[g] = xs[g] + trap(E_eval, xa_eval)/phi_g;
  }

  std::cout << "Beginning thermal group constant calculations" << std::endl;
  { // thermal group (group G+1)

    group_min = Eg[G+1]; // lower bound of group
    group_max = Eg[G]; // upper bound of group
    r = (group_max - group_min)/ng; // common difference between evaluation points
    
    // evaluation points for integral
    E_eval[0] = group_min;
    for (int i = 1; i < ng; ++i) {
      E_eval[i] = E_eval[i-1] + r;
    }
    E_eval[ng] = group_max;

    std::cout << E_eval[ng] << " should equal " << E_eval[ng-1] + r << std::endl;

    // evaluate cross sections and fluxes
    for (int i = 0; i < E_eval.size(); ++i) {

      phi_eval[i] = MaxwellBoltzmannKernel(E_eval[i]);

      dynamicNeutron->SetKineticEnergy(E_eval[i]);

      xs_eval[i] = phi_eval[i]
            *(elasticDataStore->GetCrossSection(dynamicNeutron, material) 
            + inelasticDataStore->GetCrossSection(dynamicNeutron, material));

      xa_eval[i] = phi_eval[i]*captureDataStore->GetCrossSection(dynamicNeutron,
                                  material);

    }
    phi_g = trap(E_eval, phi_eval);

    xs[G+1] = trap(E_eval, xs_eval)/phi_g;
    xt[G+1] = xs[G+1] + trap(E_eval, xa_eval)/phi_g;
  }

  std::cout << "Group constants calculated" << std::endl;

  std::string filename = output_file_base + "_" + material_name + "_xs.dat";

  std::cout << "Writing cross section data to " << filename << std::endl;

  std::ofstream outputStream(filename);

  for (int g = 0; g < G+2; ++g) {
    outputStream << g << " " << xs[g] << " " << xt[g] << "\n";
  }

  outputStream.close();

  std::cout << "Done" << std::endl;

  G4cout.rdbuf(G4cout_oldbuf);

  return 0;
}
