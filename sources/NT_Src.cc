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
  // source_weights:
  // the spectrum data in CDMSsources/spectra/neutron have 
  // units of neutrons per second per cm^3 per ppb contamination
  // 
  // we need our group source constants to have units of neutrons 
  // per second per cm^3
  // --> need to multiply by the concentration (in ppb) of each isotope
  //     (this will be stored in source_weights)
  //     (data from Table 6 of https://confluence.slac.stanford.edu/pages/viewpage.action?pageId=383932067)
  // --> also need to divide out units of energy in source spectra by dividing
  //     by bin width of spectrum
  if (material_name == "Norite" || material_name == "G4_Si" || material_name == "G4_Pb" || material_name == "G4_POLYETHYLENE" || material_name == "BoratedPoly" || material_name == "G4_Cu") {
    source_files.push_back(source_folder + "norite_2013_U_1ppb.dat");
    source_weights.push_back(1.095e3);

    source_files.push_back(source_folder + "norite_2013_Th_1ppb.dat");
    source_weights.push_back(5.715e3);

  } else {
    throw (std::runtime_error("Error in NT_Src: material not in candidate list."
    "\n    Material must be one of: Norite, "));
  }

  // double sweight_total = std::accumulate(source_weights.begin(), source_weights.end(), 0.);


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
  //Eg[G+1] = std::numeric_limits<G4double>::epsilon()*Eg[G];
  Eg[G+1] = 1e-6*Eg[G];
  

  // vector of group sources (Sg[g] refers to group g)
  doubles Sg(G+2);
  // initialize to zero
  for (int g = 0; g < G+2; ++g) {
    Sg[g] = 0;
  }


  double gmin, gmax;
  double E1, E2, dE1, dE2, s1, s2, stot;

  bool loop_again;

  std::ifstream sourcefile;
  

  for (int k = 0; k < source_files.size(); ++k) {
    // k - pull source data from file k in file list

    std::cout << "Reading spectrum from " << source_files[k] << std::endl;

    sourcefile = std::ifstream(source_files[k]);

    sourcefile >> E1 >> s1 >> E2 >> s2;

    // bin widths - need to divide "s" data by bins (assume (E2,s2) data corresponds to area of source spectrum between E1 and E2)
    dE1 = E1; // assume first "bin" in sources extends down to zero energy
    dE2 = E2 - E1;


    // if E2 < Eg[G+1], ignore the (E1,s1) ordered pair and move on
    while (E2 < Eg[G+1]) {
      E1 = E2;
      s1 = s2;
      dE1 = dE2;
      if (!sourcefile.eof()) {
        sourcefile >> E2 >> s2;
        dE2 = E2 - E1;
      } else {
        throw (std::runtime_error("The file " + source_files[k] + " does not appear to contain data in the region of interest"));
      }
    }

    for (int g = G+1; g > 0; --g) { // for each group

      gmin = Eg[g]; // lower bound of group
      gmax = Eg[g-1]; // upper bound of group

      stot = 0.;
      loop_again = true;

      do {
        if (gmin < E1) {
          if (gmax < E1) {
            loop_again = false;
            break; // skip this group
          } else {
            stot += s1;           
          }
        } else {
          stot += (E2 - gmin)*s2/(E2 - E1);
        }
        if (E2 < gmax) { // read in next values
          E1 = E2;
          dE1 = dE2;
          s1 = s2;
          if (!(sourcefile >> E2 >> s2)) {
            stot += s1;
            E1 = Eg[0] + 1.;
            E2 = Eg[0] + 1.;
            loop_again = false;
          }
        } else {
          stot += (gmax - E1)*s2/(E2 - E1);
          loop_again = false;
        }
      } while (loop_again);

      if (stot > 0) {
        Sg[g] += source_weights[k]*stot;
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



