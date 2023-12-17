// calculate boundary conditions for lab simulations using slowing-down 
// equation (assuming limiting behavior of flux out in surrounding rock is
// that of infinite-slab flux)
// 
// generates executable NT_BC
//
// V1 - 
//
// uses group cross section data in ../cross_sections/data/V1/data_<material>_<ngroups>_20_xs.dat
//
// uses group sources data in ../sources/data/V1/data_<material>_<ngroups>_Sg.dat
// 
// usage:
// ./NT_BC output_file_base_path path_to_ntransporter_base material [ngroups=100]


#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <cstdlib>

typedef std::vector<double> doubles;


int main(int argc, char *argv[]) {

    int G;
    std::string output_file_base, NT_path_base, material_name;

    // parse command line args
    if (argc < 4) {
        throw std::runtime_error("Error in NT_BC: Not enough arguments."
        "\n     Usage: ./NT_BC output_file_base_path path_to_ntransporter_base material [ngroups=100]");
    } else {
        try {
            output_file_base = argv[1];
            NT_path_base = argv[2];
            material_name = argv[3];
            G = 100; // number of fast groups
            if (argc > 4) {
                G = std::stoi(argv[4]);
            }
        } catch (std::invalid_argument) {
            throw std::runtime_error("Error in NT_BC: Invalid arguments. ngroups must be numeric."
            "\n     Usage: ./NT_BC output_file_base_path path_to_ntransporter_base material [ngroups=100]");
        }
    }

    doubles Eg(G+2), xs(G+2), xt(G+2);

    std::string xs_filename = NT_path_base + "/cross_sections/data/V1/data_" 
                        + material_name + "_" 
                        + std::to_string(G) 
                        + "_20_xs.dat";

    std::ifstream xs_stream(xs_filename);

    std::cout << "Reading cross section data from " << xs_filename << std::endl;

    int g2; // dummy variable

    for (int g = 0; g < G+2; ++g) {
        if (!(xs_stream >> g2 >> Eg[g] >> xs[g] >> xt[g])) {
            throw std::runtime_error("Error in NT_BC: reading variables for g = " + std::to_string(g) + " failed in " + xs_filename);
        }
        if (g2 != g) {
            throw std::runtime_error("Error in NT_BC: line indices in " + xs_filename + " don't match, g = " + std::to_string(g));
        }
    }

    xs_stream.close();

    std::cout << "Done reading cross section data" << std::endl;

    std::string Sg_filename = NT_path_base + "/sources/data/V1/data_" 
                        + material_name + "_" 
                        + std::to_string(G) 
                        + "_Sg.dat";

    std::ifstream Sg_stream(Sg_filename);

    std::cout << "Reading group source data from " << Sg_filename << std::endl;

    doubles Sg(G+2);

    double Eg2; // dummy variable

    // tolerance for relative error in group boundaries
    const double energy_tol = 1e-13;

    for (int g = 0; g < G+2; ++g) {
        if (!(Sg_stream >> g2 >> Eg2 >> Sg[g])) {
            throw std::runtime_error("Error in NT_BC: reading variables for g = " + std::to_string(g) + " failed in " + Sg_filename);
        }
        if (g != g2) {
            throw std::runtime_error("Error in NT_BC: line indices in " + Sg_filename + " don't match, g = " + std::to_string(g));
        }
        if (std::abs(Eg2 - Eg[g])/Eg2 > energy_tol) {
            throw std::runtime_error("Error in NT_BC: group boundary mismatch between " + xs_filename + " and " + Sg_filename + ", g = " + std::to_string(g));
        }
    }

    Sg_stream.close();

    std::cout << "Done reading group source data" << std::endl;


    std::string dx_filename = NT_path_base + "/cross_sections/data/V2/data_" 
                        + material_name + "_" 
                        + std::to_string(G) 
                        + "_dx.dat";

    std::cout << "Reading differential cross section data from " << dx_filename << std::endl;

    std::ifstream dx_stream(dx_filename);

    // variables used in loop
    int g1, gp; // counters
    double emission_density; // incoming "emission" rate
    double sig; // differential group scattering cross section
    double sigself; // self-scattering cross section
    std::stringstream line; // line stringstream
    std::string theLine;

    doubles phi(G+2); // group fluxes
    phi[0] = 0.;

    // calculate fluxes
    for (int g = 1; g < G+2; ++g) {
        
        if (!dx_stream.eof()) {
            getline(dx_stream, theLine);
        } else {
            throw std::runtime_error("Error in NT_BC: no more lines in dx file " + dx_filename);
        }

        //std::cout << "line contents pre-insertion: " << line.str() << std::endl;

        line.str(theLine);
        line.clear(); // clear error state
        
        //std::cout << "line contents post-insertion: " << line.str() << std::endl;

        //std::cout << "g = " << g << " : theLine = " << theLine << std::endl;

        line >> g1;

        //std::cout << "g1 = " << g1 << std::endl;

        if (g1 != g) {
            throw std::runtime_error("Error in NT_BC: Index mismatch at " + std::to_string(g) + " in dx file " + dx_filename);
        }

        if (!(line >> sigself)) {
            throw std::runtime_error("Error in NT_BC: no self-scattering cross section in " + dx_filename);
        } 

        //std::cout << "sigself: " << sigself << std::endl;

        emission_density = Sg[g];

        gp = g-1;

        while (line >> sig) {
            if (g <= 0) {
                throw std::runtime_error("Error in NT_BC: Too many entries in line " + std::to_string(g) + " of " + dx_filename);
            }
            emission_density += sig*phi[gp];
            --gp;
        }

        //phi[g] = (xs[g-1]*phi[g-1] + Sg[g])/xt[g];
        phi[g] = emission_density/(xt[g] - sigself);
    }
    
    dx_stream.close();

    std::string output_filename = output_file_base + "_"
                        + material_name + "_"
                        + std::to_string(G) 
                        + "_BC_V2.dat";

    std::cout << "Writing flux data to " << output_filename << std::endl;

    std::ofstream output_stream(output_filename);

    output_stream << std::setprecision(17);

    for (int g = 0; g < G+2; ++g) {
        output_stream << g << " " << Eg[g] << " " << phi[g] << "\n";
    }

    output_stream.close();

    std::cout << "Done" << std::endl;

    return 0;
}




