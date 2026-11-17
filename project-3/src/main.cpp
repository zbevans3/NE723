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

void reset_solution(coreData& core,
        quadratureData& quadrature,
        meshData& mesh,
        solverData& solver,
        QDData& QD,
        QDMGData& QDMG,
        QD1GData& QD1G);

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
    QDData QD;
    QDMGData QDMG;
    QD1GData QD1G;
    quadratureData quadrature;
    std::vector<materialData> materials;

    try {
        //reading input file
        std::cout << "Reading input file\n";

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
    std::cout << "Building mesh\n";    
    Mesh buildmesh;
    buildmesh.build(core, mesh);
    if (core.search != "critical_size") out.write_mesh(outfile, mesh, materials);

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
    reset_solution(core, quadrature, mesh, solver, QD, QDMG, QD1G);

    // start solve
    if (core.search == "fixed_source"){
       // fixed source calculation
       std::cout << "Performing fixed source calculation\n";

       solver.lambda = 1.0;
       double invlambda = solver.lambda;

       SourceIteration inner;
       inner.solve(core, 
                   quadrature,
                   materials,
                   mesh,
                   solver,
                   QD,
                   QDMG,
                   QD1G,
                   invlambda);

       out.write_inner_iteration(outfile, core, solver);

    } else if (core.search == "eigenvalue"){
       // eigenvalue calculation
       std::cout << "Performing eigenvalue calculation\n";
       bool foundfis = false;
       for (int i = 0 ; i < mesh.ncells ; i++){
          if (materials[mesh.mapmat[i]].xsfiss.back() > 0.0){
            // check if most thermal group has fission
            foundfis = true;
            break;
          }
       }
 
       if (foundfis == false){
        std::cerr << "ERROR: No fissionable material found in mesh. Exiting...\n";
        return 1;
       }

       solver.lambda = 1.0;
       OuterIteration outer;
       outer.solve(core,
                   quadrature,
                   materials,
                   mesh,
                   solver,
                   QD,
                   QDMG,
                   QD1G);

       out.write_outer_iteration(outfile, core, solver);

    } else if (core.search == "critical_size"){
       std::cout << "Performing critical size search\n";
       if (core.nregions > 1){
          std::cerr << "ERROR: Only one region can be used in critical size search. Exiting...\n";
          return 1;
       }

       double hcrit = core.xedge[0];  // initial guess comes from input
       int lcrit = 0;
       double kdiff = 0.0;
       double hprev = hcrit;
       double kdiffprev = kdiff;
       double hnew = hcrit;

       core.itdebug = false;    // suppress inner iteration priting 
               
       // perform initial solve
       solver.lambda = 1.0;
       OuterIteration outer;
       outer.solve(core,
                  quadrature,
                  materials,
                  mesh,
                  solver,
                  QD,
                  QDMG,
                  QD1G);
                  
       kdiff = solver.lambda - core.kcrit;   

       std::cout << std::right
                 << std::setw(12) << "Iteration"
                 << std::setw(16) << "h [cm]"
                 << std::setw(16) << "keff"
                 << std::setw(16) << "kcrit"
                 << std::setw(16) << "difference"
                 << '\n';   
       std::cout << std::right
                 << std::setw(12) << lcrit
                 << std::fixed << std::setprecision(8)
                 << std::setw(16) << hcrit
                 << std::setw(16) << solver.lambda
                 << std::setw(16) << core.kcrit
                 << std::scientific << std::setprecision(4)
                 << std::setw(16) << kdiff
                 << '\n';   
       outfile << " Critical size search results:\n";
       outfile << std::right
                 << std::setw(12) << "Iteration"
                 << std::setw(16) << "h [cm]"
                 << std::setw(16) << "keff"
                 << std::setw(16) << "kcrit"
                 << std::setw(16) << "difference"
                 << '\n';
       outfile << "   " << std::string(73, '-') << '\n';   
       outfile << std::right
                 << std::setw(12) << lcrit
                 << std::fixed << std::setprecision(8)
                 << std::setw(16) << hcrit
                 << std::setw(16) << solver.lambda
                 << std::setw(16) << core.kcrit
                 << std::scientific << std::setprecision(4)
                 << std::setw(16) << kdiff
                 << '\n';           

       // larger if subcritical, smaller if supercritical.
       if (std::fabs(kdiff) > core.epsk) {
           hcrit *= (kdiff < 0.0) ? 1.1 : 0.9;
       }
       
       while (true) {
       
           ++lcrit;
           core.xedge[0] = hcrit;
       
           // Build and solve at the new size.
           buildmesh.build(core, mesh);
           reset_solution(core, quadrature, mesh, solver, QD, QDMG, QD1G);
       
           solver.lambda = 1.0;
           outer.solve(core,
                       quadrature,
                       materials,
                       mesh,
                       solver,
                       QD,
                       QDMG,
                       QD1G);
       
           kdiff = solver.lambda - core.kcrit;
       
           std::cout << std::right
                     << std::setw(12) << lcrit
                     << std::fixed << std::setprecision(8)
                     << std::setw(16) << hcrit
                     << std::setw(16) << solver.lambda
                     << std::setw(16) << core.kcrit
                     << std::scientific << std::setprecision(4)
                     << std::setw(16) << kdiff
                     << '\n';   
           outfile << std::right
                     << std::setw(12) << lcrit
                     << std::fixed << std::setprecision(8)
                     << std::setw(16) << hcrit
                     << std::setw(16) << solver.lambda
                     << std::setw(16) << core.kcrit
                     << std::scientific << std::setprecision(4)
                     << std::setw(16) << kdiff
                     << '\n';                      
       
           if (std::fabs(kdiff) <= core.epsk) {
               break;
           }
       
           // Secant estimate of the critical size.
           hnew = hcrit - kdiff * (hcrit - hprev) / (kdiff - kdiffprev);
       
           // Keep the trial size positive.
           if (hnew <= 0.0) {
               hnew = 0.5 * hcrit;
           }
       
           hprev = hcrit;
           kdiffprev = kdiff;
           hcrit = hnew;
       }
       outfile << "   hcrit: " << std::fixed << std::setprecision(8)
                             << std::setw(16) << hcrit << " cm\n";
       outfile << "\n";

       // write final mesh
       outfile << " Final mesh:\n";
       out.write_mesh(outfile, mesh, materials);

    } else {
        std::cerr << "Unknown search type. Exiting...\n";
        return 1;
    }

    out.write_angular_flux(outfile, mesh, quadrature, solver);
    out.write_scalar_flux(outfile, mesh, solver);
    out.write_current(outfile, mesh, solver);
    out.write_residuals(outfile, mesh, solver);    

    if (core.acceleration != "none"){
        if (core.acceleration == "QD"){
            // write LO solution data
            out.write_QD_data(outfile, mesh, QD);
            out.write_scalar_flux_LO(outfile, mesh, QD);
            out.write_current_LO(outfile, mesh, QD);
        }
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    time = duration.count() / 1e6;
    (exitF == true) ? std::cout << "Execution finished successfully!\n" : std::cout << "Execution failed.\n";
    std::cout<<"Execution time (sec):  " << std::fixed << std::setprecision(4) << time << "\n";

    outfile.close();
    return 0;
}

void reset_solution(coreData& core,
        quadratureData& quadrature,
        meshData& mesh,
        solverData& solver,
        QDData& QD,
        QDMGData& QDMG,
        QD1GData& QD1G){

    int ngroups = core.ngroups;
    double leftval  = (mesh.mapbc.front() == 3) ? core.psileft  : 0.0;
    double rightval = (mesh.mapbc.back()  == 3) ? core.psiright : 0.0;
    std::vector<std::vector<double>> psi_edge(quadrature.order, std::vector<double>(mesh.nedges, 1.0));
    std::vector<std::vector<double>> psi_cell(quadrature.order, std::vector<double>(mesh.ncells, 1.0));
    solver.psie.assign(ngroups, psi_edge);
    solver.psibar.assign(ngroups, psi_cell);
    solver.psihat.assign(ngroups, psi_cell);
    solver.psilbc.assign(ngroups, std::vector<double>(quadrature.order, leftval));
    solver.psirbc.assign(ngroups, std::vector<double>(quadrature.order, rightval));
    solver.phie.assign(ngroups, std::vector<double>(mesh.nedges, 1.0));
    solver.phibar.assign(ngroups, std::vector<double>(mesh.ncells, 1.0));
    solver.phihat.assign(ngroups, std::vector<double>(mesh.ncells, 1.0));
    solver.je.assign(ngroups, std::vector<double>(mesh.nedges, 0.0));
    solver.jbar.assign(ngroups, std::vector<double>(mesh.ncells, 0.0));
    solver.jhat.assign(ngroups, std::vector<double>(mesh.ncells, 0.0));
    solver.source.assign(ngroups, std::vector<double>(mesh.ncells, 0.0));
    solver.source1.assign(ngroups, std::vector<double>(mesh.ncells, 0.0));
    solver.sresidual.assign(ngroups, std::vector<double>(mesh.ncells, 0.0));
    solver.sresidual1.assign(ngroups, std::vector<double>(mesh.ncells, 0.0));
    solver.iresidual.assign(ngroups, std::vector<double>(mesh.ncells, 0.0));
    solver.iresidual1.assign(ngroups, std::vector<double>(mesh.ncells, 0.0));
    solver.keff.clear();    
    solver.ispecrad.clear();
    solver.ospecrad.clear();
    solver.ilinfphi.clear();
    solver.olinfphi.clear();
    solver.linffsrc.clear();
    if (core.acceleration == "QD"){
        // build QD data
        QD.CbfL = -0.5;
        QD.CbfR = 0.5;
        QD.Edde.assign(mesh.nedges, 1.0/3.0);
        QD.Eddbar.assign(mesh.ncells, 1.0/3.0);
        QD.phibarLO.assign(mesh.ncells, 0.0);
        QD.phieLO.assign(mesh.nedges, 0.0);
        QD.jeLO.assign(mesh.nedges, 0.0);
        QD.mLO.assign(mesh.ncells, 0.0);
        QD.mprimeLO.assign(mesh.ncells, 0.0);
    } else if (core.acceleration == "MGQD"){
        // build ML QD data
        QDMG.Edde.assign(ngroups, std::vector<double>(mesh.nedges, 1.0/3.0));
        QDMG.Eddbar.assign(ngroups, std::vector<double>(mesh.ncells, 1.0/3.0));
        QDMG.phie.assign(ngroups, std::vector<double>(mesh.nedges, 1.0));
        QDMG.phibar.assign(ngroups, std::vector<double>(mesh.ncells, 1.0));
        QDMG.je.assign(ngroups, std::vector<double>(mesh.nedges, 1.0));
        QDMG.CbfL.assign(ngroups, 0.0);
        QDMG.CbfR.assign(ngroups, 0.0);

        QD1G.phie.assign(mesh.nedges, 1.0);
        QD1G.phibar.assign(mesh.ncells, 1.0);
        QD1G.je.assign(mesh.nedges, 1.0);
        QD1G.xsabs.assign(mesh.ncells, 0.0);
        QD1G.xsnusigf.assign(mesh.ncells, 0.0);
        QD1G.Hm.assign(mesh.nedges, 1.0/3.0);
        QD1G.Hp.assign(mesh.nedges, 1.0/3.0);
        QD1G.Am.assign(mesh.nedges, 1.0/3.0);
        QD1G.Ap.assign(mesh.nedges, 1.0/3.0);
        QD1G.CbfL = -0.5;
        QD1G.CbfR = 0.5;

    }

}