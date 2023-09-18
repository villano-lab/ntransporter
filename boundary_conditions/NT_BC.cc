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


    // calculate fluxes
    doubles phi(G+2);
    phi[0] = 0.;

    for (int g = 1; g < G+2; ++g) {
        phi[g] = (xs[g-1]*phi[g-1] + Sg[g])/xt[g];
    }

    std::string output_filename = output_file_base + "_"
                        + material_name + "_"
                        + std::to_string(G) 
                        + "_BC_V1.dat";

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




