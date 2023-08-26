// calculate group source terms Sg for radiogenic sources in specified material
// 
// generates executable NT_Src
//
// usage:
// ./NT_Src output_file_base_path material path_to_supersim [ngroups=100]

#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <cmath>
#include <limits>

#include "NTUtilities.hh"

typedef std::vector<double> doubles;


int main(int argc, char *argv[]) {

  int G;
  std::string output_file_base, material_name, path_to_supersim;

  // parse command line args
  if (argc < 4)  {
    throw std::runtime_error("Error in NT_Src: Not enough arguments."
         "\n    Usage: ./NT_Src output_file_base_path material path_to_supersim " 
         "[ngroups=100] [points_per_group=10]");
  } else {
    try {
      output_file_base = argv[1];
      material_name = argv[2];
      path_to_supersim = argv[3];
      G = 100; // number of groups
      if (argc > 4) {
        G = std::stoi(argv[4]);
      }
    } catch (std::invalid_argument) {
      throw std::runtime_error("Error in NT_Src: Invalid arguments. ngroups " 
        "must be numeric."
        "\n    Usage: ./NT_Src output_file_base_path material path_to_supersim "
        "[ngroups=100] [points_per_group=10]");
    }
  }
  
  // folder containing neutron spectrum data in SuperSim
  std::string source_folder = path_to_supersim + "/CDMSsources/spectra/neutron/";


  std::vector<std::string> source_files;
  doubles source_weights;
  if (material_name == "Norite") {
    source_files.push_back(source_folder + "norite_2013_U_1ppb.dat");
    source_weights.push_back(1.);

    source_files.push_back(source_folder + "norite_2013_Th_1ppb.dat");
    source_weights.push_back(1.);

  } else {
    throw (std::runtime_error("Error in NT_Src: material not in candidate list."
    "\n    Material must be one of: Norite, "))
  }

  // alpha (common ratio of group boundaries)
  double alpha = std::pow(Emin/Emax, 1./G);

  // array of group boundaries (one thermal group)
  doubles Eg(G+2);
  Eg[0] = Emax;

  // calculate group boundaries
  for (int g = 1; g < G+1; ++g) {
    Eg[g] = Eg[g-1]*alpha;
  }
  // lower bound (basically zero)
  Eg[G+1] = std::numeric_limits<G4double>::epsilon()*Eg[G];
  

  // vector of group sources (first element always zero; Sg[g] refers to group g)
  doubles Sg(G+2);
  Sg[0] = 0.;


  
  double group_min, group_max, r, phi_g;
  
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
            *getSourceStrength(E_eval[i]);

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

    //std::cout << E_eval[ng-1] << " should equal " << E_eval[ng] - r << std::endl;

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

  return 0;
}



