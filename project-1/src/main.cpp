#include "data.h"
#include "write_output.h"
#include "read_input.h"
#include "build_mesh.h"
#include "source_iteration.h"
#include "outer_iteration.h"

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
    bool exitF = true;

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
    out.write_mesh(outfile, core, mesh, materials);

    // build reflective ordinate map
    if (core.bcleft == "reflective" || core.bcright == "reflective") {
        quadrature.reflmap.assign(quadrature.order, -1);
        for (int u1 = 0; u1 < quadrature.order; ++u1) {
            for (int u2 = 0; u2 < quadrature.order; ++u2) {
                if (std::fabs(quadrature.angle[u1] +
                              quadrature.angle[u2]) < 1.0e-10) {
                    quadrature.reflmap[u1] = u2;
                    break;
                }
            }
        }
    }

    // allocate solver arrays
    solver.psie.assign(quadrature.order, std::vector<double>(mesh.nedges, 1.0));
    solver.psibar.assign(quadrature.order, std::vector<double>(mesh.ncells, 1.0));
    solver.psihat.assign(quadrature.order, std::vector<double>(mesh.ncells, 1.0));
    solver.phie.assign(mesh.nedges, 1.0);
    solver.phibar.assign(mesh.ncells, 1.0);
    solver.phihat.assign(mesh.ncells, 1.0);
    solver.je.assign(mesh.nedges, 0.0);
    solver.jbar.assign(mesh.ncells, 0.0);
    solver.jhat.assign(mesh.ncells, 0.0);
    solver.source.assign(mesh.ncells, 0.0);
    solver.source1.assign(mesh.ncells, 0.0);
    solver.sresidual.assign(mesh.ncells, 0.0);
    solver.sresidual1.assign(mesh.ncells, 0.0);
    solver.iresidual.assign(mesh.ncells, 0.0);
    solver.iresidual1.assign(mesh.ncells, 0.0);

    // start solve
    if (core.search == "fixed_source"){
       // fixed source calculation
       std::cout << "performing fixed source calculation\n";
       std::cout << "starting source iteration\n";

       solver.lambda = 1.0;
       double invlambda = solver.lambda;

       SourceIteration inner;
       inner.solve(core, quadrature, materials, mesh, solver, invlambda);

    } else if (core.search == "eigenvalue"){
       // eigenvalue calculation
       std::cout << "performing eigenvalue calculation\n";
       bool foundfis = false;
       for (int i = 0 ; i < mesh.ncells ; i++){
          if (materials[mesh.mapmat[i]].xsfiss > 0.0){
            foundfis = true;
            break;
          }
       }
 
       if (foundfis == false){
        exitF = false;
        std::cerr << "ERROR: No fissionable material found in mesh. Exiting.\n";
       }

       solver.lambda = 1.0;
       OuterIteration outer;
       outer.solve(core, quadrature, materials, mesh, solver);

    } else if (core.search == "critical_size"){
       std::cout << "performing critical size search\n";
       if (core.nregions > 1){
          exitF = false;            
          std::cerr << "ERROR: Only one region can be used in critical size search. Exiting.\n";
       }

    } else {
        std::cerr << "Unknown type. Exiting...\n";
        return 1;
    }

    out.write_angular_flux(outfile, mesh, quadrature, solver);
    out.write_scalar_flux(outfile, mesh, quadrature, solver);
    out.write_current(outfile, mesh, quadrature, solver);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    time = duration.count() / 1e6;
    (exitF == true) ? std::cout << "Execution finished successfully!\n" : std::cout << "Execution failed.\n";
    std::cout<<"Execution time (sec):  " << std::fixed << std::setprecision(4) << time << "\n";

    outfile.close();
    return 0;
}

