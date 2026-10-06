#include "data.h"
#include "write_output.h"
#include "read_input.h"
#include "build_mesh.h"

#include <exception>
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <chrono>

int main(int argc, char* argv[])
{
    double time;
    auto start = std::chrono::high_resolution_clock::now();

    std::string inpfile;
    inpfile=argv[1];

    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input file>\n";
        return 1;
    }

    // opening output file and write header
    Output out;
    std::string outputFile = inpfile.substr(0, inpfile.length()-4);
    outputFile += ".out";
    std::ofstream outfile(outputFile);
    out.write_header(outfile);

    coreData core;
    meshData mesh;
    solverData solver;
    quadratureData quadrature;
    std::vector<materialData> materials;

    try {
        //reading input file
        std::cout << "reading input file\n";

        // read input
        Input reader;
        reader.read(inpfile, core, quadrature, materials);

        // echo input
        out.echo_input(outfile, core, quadrature, materials);

    } catch (const std::exception& error) {
        std::cerr << "Input error: " << error.what() << '\n';
        return 1;
    }

    // build mesh
    std::cout << "building mesh\n";    
    Mesh buildmesh;
    buildmesh.build(core, mesh);
    out.write_mesh(outfile, core, mesh);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    time = duration.count() / 1e6;
    std::cout << "Execution finished successfully!\n";
    std::cout<<"Execution time (sec):  " << std::fixed << std::setprecision(4) << time << "\n";

    outfile.close();
    return 0;
}

