// calculate boundary conditions for lab simulations using slowing-down 
// equation (assuming limiting behavior of flux out in surrounding rock is
// that of infinite-slab flux)
// 
// generates executable NT_BC
// 
// usage:
// ./NT_BC output_file_base_path material [ngroups=100]


#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <cmath>

typedef std::vector<double> doubles;


int main(int argc, char *argv[]) {

    int G;
    std::string output_file_base, material_name;

    // parse command line args
    if (argc < 3) {
        throw std::runtime_error("Error in NT_BC: Not enough arguments."
        "\n     Usage: ./NT_BC output_file_base_path material [ngroups=100]");
    } else {
        try {
            output_file_base = argv[1];
            material_name = argv[2];
            G = 100; // number of fast groups
            if (argc > 3) {
                G = std::stoi(argv[3]);
            }
        } catch (std::invalid_argument) {
            throw std::runtime_error("Error in NT_BC: Invalid arguments. ngroups must be numeric."
            "\n     Usage: ./NT_BC output_file_base_path material [ngroups=100]");
        }
    }

    doubles Eg(G+2), xs(G+2), xt(G+2);

    std::string xs_filename = "../cross_sections/data/V1/data_" 
                          + material_name + "_" 
                          + std::to_string(G) 
                          + "_20_xs.dat";

    std::ifstream xs_stream(xs_filename);

    std::cout << "Reading cross section data from " << xs_filename << std::endl;

    for (int g = 0; g < G+2; ++g) {
        xs_stream >> Eg[g] >> xs[g] >> xt[g];
    }

    xs_stream.close();


    return 0;
}




