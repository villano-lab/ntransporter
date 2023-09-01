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
#include <numeric>
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
         "[ngroups=100]");
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
        "[ngroups=100]");
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
    "\n    Material must be one of: Norite, "));
  }

  double sweight_total = std::accumulate(source_weights.begin(), source_weights.end(), 0.);


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
  Eg[G+1] = std::numeric_limits<double>::epsilon()*Eg[G];


  // vectors of energies, source rates, and fluxes for which there is data in the .dat files
  doubles E_eval{0.};
  doubles S_eval{0.};
  doubles phi_eval{0.};
  

  // vector of group sources (Sg[g] refers to group g)
  doubles Sg(G+2);
  // initialize to zero
  for (int g = 0; g < G+2; ++g) {
    Sg[g] = 0;
  }


  double gmin, gmax, r, phi_g;
  double E1, E2, s1, s2;

  bool calc_this;

  std::ifstream sourcefile;


  double (*phi_func)(double);
  

  for (int k = 0; k < source_files.size(); ++k) {
    // k - pull source data from file k in file list

    std::cout << "Reading spectrum from " << source_files[k] << std::endl;

    sourcefile = std::ifstream(source_files[k]);

    sourcefile >> E1 >> s1 >> E2 >> s2;


    // if E2 < Eg[G+1], ignore the (E1,s1) ordered pair and move on
    while (E2 < Eg[G+1]) {
      E1 = E2;
      s1 = s2;
      if (!sourcefile.eof()) {
        sourcefile >> E2 >> s2;
      } else {
        throw std::exception("The file " source_files[k] " does not appear to contain data in the region of interest");
      }
    }

    for (int g = G+1; g > 0; --g) { // for each group
      if (g%10 == 1) {
        std::cout << (G+2-g) << " / " << G+1 << std::endl;
      }
      if (g == G+1) { // thermal group
        phi_func = MaxwellBoltzmannKernel;
      } else { // fast group
        phi_func = fast_kernel;
      }


      gmin = Eg[g]; // lower bound of group
      gmax = Eg[g-1]; // upper bound of group

      E_eval.clear();
      phi_eval.clear();
      S_eval.clear();
      calc_this = true;

      do {
        if (gmin < E1) {
          if (gmax < E1) {
            calc_this = false; 
            break; // skip this group
          } else {
            E_eval.push_back(E1);
            phi_eval.push_back(phi_func(E1));
            S_eval.push_back(s1*phi_eval.back());            
          }
        } else {
          E_eval.push_back(gmin);
          phi_eval.push_back(phi_func(gmin));
          S_eval.push_back(phi_eval.back()*interp(gmin, E1, s1, E2, s2));
        }
        if (E2 < gmax) { // read in next values
          E1 = E2;
          s1 = s2;
          if (!sourcefile.eof()) {
            sourcefile >> E2 >> s2;
          } else {
            E2 = gmax;
          }
        } else {
          E_eval.push_back(gmax);
          phi_eval.push_back(phi_func(gmax));
          S_eval.push_back(phi_eval.back()*interp(gmax, E1, s1, E2, s2));
        }
      } while (E2 < gmax);

      if (calc_this) {
        phi_g = trap(E_eval, phi_eval);
        Sg[g] += (source_weights[k]/sweight_total)*trap(E_eval, S_eval)/phi_g;
      }
    }
  }

  std::cout << "All group constants calculated" << std::endl;

  std::string filename = output_file_base + "_" 
                       + material_name + "_" 
                       + std::to_string(G)
                       + "_Sg.dat";

  std::cout << "Writing group source data to " << filename << std::endl;

  std::ofstream outputStream(filename);

  outputStream << std::setprecision(17);

  for (int g = 0; g < G+2; ++g) {
    outputStream << g << " " << Eg[g] << " "<< Sg[g] << "\n";
  }

  outputStream.close();

  std::cout << "Done" << std::endl;

  return 0;
}



