#include "data.h"
#include "read_input.h"

#include <exception>
#include <iostream>

int main(int argc, char* argv[])
{
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input file>\n";
        return 1;
    }

    try {
        coreData core;
        solverData solver;
        quadratureData quadrature;
        std::vector<materialData> materials;

        Input reader;
        reader.read(argv[1], core, materials);

        // core and materials now contain the parsed input.
        // Allocate solver arrays and construct quadrature here.

        std::cout << "Materials: " << core.nmat << '\n';
        std::cout << "Material 1 total XS: "
                  << materials[0].xstot << '\n';

    } catch (const std::exception& error) {
        std::cerr << "Input error: " << error.what() << '\n';
        return 1;
    }

    return 0;
}

