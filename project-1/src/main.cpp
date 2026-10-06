#include "data.h"
#include "read_input.h"
#include "build_mesh.h"

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
        meshData mesh;
        solverData solver;
        quadratureData quadrature;
        std::vector<materialData> materials;

        // read input
        Input reader;
        reader.read(argv[1], core, quadrature, materials);

        // build mesh
        Mesh buildmesh;
        buildmesh.build(core, mesh);

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

