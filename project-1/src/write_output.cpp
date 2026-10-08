#include "write_output.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <numeric>

void Output::write_header(std::ofstream &outfile)
{
    outfile << "===================================================================\n";
    outfile << "=                            NE723 Project 1                      =\n";
    outfile << "=                            Zachary Bevans                       =\n";
    outfile << "=                              10/13/2026                         =\n";
    outfile << "===================================================================\n";
    outfile << "Code to solve the steady state one group neutron transport equation\n";
    outfile << " in 1D slab geometry the discrete ordinates method.\n";
    outfile << "===================================================================\n\n";
}

void Output::echo_input(
    std::ofstream &outfile,
    coreData& core,
    quadratureData& quadrature,
    std::vector<materialData>& materials
){
   
    outfile << " Input values:\n";
    outfile << "  [CORE]\n";
    outfile << std::left;
    outfile << "    " << std::setw(12) << "search"   << core.search   << "\n";
    outfile << "    " << std::setw(12) << "itdebug"  << core.itdebug  << "\n";
    outfile << "    " << std::setw(12) << "epsk"     << core.epsk     << "\n";
    outfile << "    " << std::setw(12) << "epsflx"   << core.epsflx   << "\n";
    outfile << "    " << std::setw(12) << "maxin"    << core.maxin    << "\n";
    outfile << "    " << std::setw(12) << "maxout"   << core.maxout   << "\n";
    outfile << "    " << std::setw(12) << "nregions" << core.nregions << "\n";
    outfile << "    " << std::setw(12) << "nmat"     << core.nmat     << "\n";
    outfile << "    " << std::setw(12) << "bcleft"   << std::setw(12) << core.bcleft   
            << std::setw(6) << std::fixed << std::setprecision(3) << core.psileft  << "\n";
    outfile << "    " << std::setw(12) << "bcright"   << std::setw(12) <<core.bcright   
            << std::setw(5) << std::fixed << std::setprecision(3) << core.psiright  << "\n";
    outfile << std::right;
    outfile << "    nx ";
    for (int i = 0 ; i < core.nregions ; i++){
        outfile << std::setw(6) << core.nx[i];
    }
    outfile << "\n";
    outfile << "    xedge ";
    for (int i = 0 ; i < core.nregions ; i++){
        outfile << std::setw(6) << std::fixed << std::setprecision(3) << core.xedge[i];
    }
    outfile << "\n\n";

    outfile << "  [MATERIAL]\n";
    for (int i = 0 ; i < core.nmat ; i++){
        outfile << "    material " << std::setw(6) << materials[i].id << std::setw(10) 
                << std::fixed << std::setprecision(4) << materials[i].xstot << std::setw(10)
                << std::fixed << std::setprecision(4) << materials[i].xsscat << std::setw(10)
                << std::fixed << std::setprecision(4) << materials[i].xsfiss << std::setw(10)
                << std::fixed << std::setprecision(4) << materials[i].nu << std::setw(10)
                << std::fixed << std::setprecision(4) << materials[i].q << "\n";
    }
    outfile << "\n";

    outfile << "  [QUADRATURE]\n";
    outfile << "    order " << std::setw(6) << std::right << quadrature.order << "\n";
    outfile << "    angle ";
    for (int i = 0 ; i < quadrature.order ; i++){
        outfile << std::setw(10) << std::right << std::fixed << std::setprecision(6) << quadrature.angle[i];
    }
    outfile << "\n";
    outfile << "    weight ";
    for (int i = 0 ; i < quadrature.order ; i++){
        outfile << std::setw(10) << std::right << std::fixed << std::setprecision(6) << quadrature.weight[i];
    }
    outfile << "\n";
    outfile << "    sum of weights "<< std::setw(10) << std::right << std::fixed << std::setprecision(6) 
            << std::accumulate(quadrature.weight.begin(), quadrature.weight.end(), 0.0) << "\n";
    double wsum = 0.0;
    for (int i = 0 ; i < quadrature.order ; i++){
        wsum += quadrature.angle[i]*quadrature.weight[i];
    }
    outfile << "    sum of weight*angle "<< std::setw(10) << std::right << std::fixed << std::setprecision(6) << wsum << "\n";    
    outfile << "\n\n";
}

void Output::write_mesh(
    std::ofstream& outfile,
    coreData& core,
    meshData& mesh,
    std::vector<materialData>& materials)
{

    outfile << " Mesh information:\n";
    outfile << "    Number of cells: " << mesh.ncells << '\n';

    outfile << std::right
            << std::setw(8)  << "Cell"
            << std::setw(10) << "Material"
            << std::setw(15) << "Left edge"
            << std::setw(15) << "Right edge"
            << std::setw(15) << "Center"
            << std::setw(15) << "Width"
            << '\n';

    outfile << "    " << std::string(74, '-') << '\n';
    outfile << std::fixed << std::setprecision(6);

    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(8)  << i + 1
                << std::setw(10) << materials[mesh.mapmat[i]].id
                << std::setw(15) << mesh.xloc[i]
                << std::setw(15) << mesh.xloc[i + 1]
                << std::setw(15) << mesh.xcenter[i]
                << std::setw(15) << mesh.hx[i]
                << '\n';
    }

    outfile << "\n Edge information:\n";
    outfile << "    Number of edges: " << mesh.nedges << "\n";
    outfile << std::setw(8)  << "Edge" << std::setw(15) << "Location" << std::setw(10) << "ibc" << '\n';
    outfile << "    " << std::string(29, '-') << '\n';
    for (int i = 0; i < mesh.nedges; ++i) {
        outfile << std::setw(8)  << i + 1 << std::setw(15) << mesh.xloc[i] << std::setw(10) << mesh.mapbc[i] << '\n';
    }
    outfile << '\n';
}

void Output::write_angular_flux(
    std::ofstream& outfile,
    meshData& mesh,
    quadratureData& quadrature,
    solverData& solver)
{
    const std::string indent = "    ";
    constexpr int label_width = 8;
    constexpr int cell_width = 8;
    const int table_width = label_width + cell_width * mesh.ncells;
    const int table_width_e = label_width + cell_width * mesh.nedges;

    outfile << " Cell average angular flux solution:\n";
    outfile << indent << std::right
            << std::setw(label_width) << "Angle";
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width, '-') << '\n';
    outfile << std::fixed << std::setprecision(3);
    for (int u = 0; u < quadrature.order; ++u) {
        outfile << indent
                << std::setw(label_width) << quadrature.angle[u];
        for (int i = 0; i < mesh.ncells; ++i) {
            outfile << std::setw(cell_width) << solver.psibar[u][i];
        }
        outfile << '\n';
    }
    outfile << '\n';

    outfile << " First Legendre spatial moment of the angular flux:\n";
    outfile << indent << std::right
            << std::setw(label_width) << "Angle";
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width, '-') << '\n';
    outfile << std::fixed << std::setprecision(3);
    for (int u = 0; u < quadrature.order; ++u) {
        outfile << indent
                << std::setw(label_width) << quadrature.angle[u];
        for (int i = 0; i < mesh.ncells; ++i) {
            outfile << std::setw(cell_width) << solver.psihat[u][i];
        }
        outfile << '\n';
    }
    outfile << '\n';

    outfile << " Cell edge angular flux solution:\n";
    outfile << indent << std::right
            << std::setw(label_width) << "Angle";
    for (int i = 0; i < mesh.nedges; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width_e, '-') << '\n';
    outfile << std::fixed << std::setprecision(3);
    for (int u = 0; u < quadrature.order; ++u) {
        outfile << indent
                << std::setw(label_width) << quadrature.angle[u];
        for (int i = 0; i < mesh.nedges; ++i) {
            outfile << std::setw(cell_width) << solver.psie[u][i];
        }
        outfile << '\n';
    }
    outfile << '\n';

}

void Output::write_scalar_flux(
    std::ofstream& outfile,
    meshData& mesh,
    quadratureData& quadrature,
    solverData& solver)
{
    const std::string indent = "    ";
    constexpr int label_width = 8;
    constexpr int cell_width = 8;
    const int table_width = cell_width * mesh.ncells;
    const int table_width_e = cell_width * mesh.nedges;

    outfile << " Cell average scalar flux solution:\n";
    outfile << indent << std::right
            << std::setw(label_width);
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width, '-') << '\n';
    outfile << indent << std::setw(label_width)
            << std::fixed << std::setprecision(3);
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << solver.phibar[i];
    }
    outfile << "\n\n";

    outfile << " First Legendre spatial moment of the scalar flux:\n";
    outfile << indent << std::right
            << std::setw(label_width);
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width, '-') << '\n';
    outfile << indent << std::setw(label_width)
            << std::fixed << std::setprecision(3);
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << solver.phihat[i];
    }
    outfile << "\n\n";    

    outfile << " Cell edge scalar flux solution:\n";
    outfile << indent << std::right
            << std::setw(label_width);
    for (int i = 0; i < mesh.nedges; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width_e, '-') << '\n';
    outfile << indent << std::setw(label_width)
            << std::fixed << std::setprecision(3);
    for (int i = 0; i < mesh.nedges; ++i) {
        outfile << std::setw(cell_width) << solver.phie[i];
    }
    outfile << "\n\n";
}

void Output::write_current(
    std::ofstream& outfile,
    meshData& mesh,
    quadratureData& quadrature,
    solverData& solver)
{
    const std::string indent = "    ";
    constexpr int label_width = 8;
    constexpr int cell_width = 8;
    const int table_width = cell_width * mesh.ncells;
    const int table_width_e = cell_width * mesh.nedges;

    outfile << " Cell average current solution:\n";
    outfile << indent << std::right
            << std::setw(label_width);
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width, '-') << '\n';
    outfile << indent << std::setw(label_width)
            << std::fixed << std::setprecision(3);
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << solver.jbar[i];
    }
    outfile << "\n\n";

    outfile << " First Legendre spatial moment of the current:\n";
    outfile << indent << std::right
            << std::setw(label_width);
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width, '-') << '\n';
    outfile << indent << std::setw(label_width)
            << std::fixed << std::setprecision(3);
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << solver.jhat[i];
    }
    outfile << "\n\n";    

    outfile << " Cell edge current solution:\n";
    outfile << indent << std::right
            << std::setw(label_width);
    for (int i = 0; i < mesh.nedges; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width_e, '-') << '\n';
    outfile << indent << std::setw(label_width)
            << std::fixed << std::setprecision(3);
    for (int i = 0; i < mesh.nedges; ++i) {
        outfile << std::setw(cell_width) << solver.je[i];
    }
    outfile << "\n\n";
}

void Output::write_residuals(
    std::ofstream& outfile,
    meshData& mesh,
    quadratureData& quadrature,
    solverData& solver)
{
    const std::string indent = "    ";
    constexpr int label_width = 14;
    constexpr int cell_width = 14;
    const int table_width = cell_width * mesh.ncells;

    outfile << " Cell average scheme residual:\n";
    outfile << indent << std::right
            << std::setw(label_width);
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width, '-') << '\n';
    outfile << indent << std::setw(label_width)
            << std::scientific << std::setprecision(4);
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << solver.sresidual[i];
    }
    outfile << "\n\n";

    outfile << " Scheme residual for first Legendre spatial moment:\n";
    outfile << indent << std::right
            << std::setw(label_width);
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width, '-') << '\n';
    outfile << indent << std::setw(label_width)
            << std::scientific << std::setprecision(4);
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << solver.sresidual1[i];
    }
    outfile << "\n\n";

    outfile << " Cell average iterative residual:\n";
    outfile << indent << std::right
            << std::setw(label_width);
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width, '-') << '\n';
    outfile << indent << std::setw(label_width)
            << std::scientific << std::setprecision(4);
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << solver.sresidual[i];
    }
    outfile << "\n\n";

    outfile << " Iterative residual for first Legendre spatial moment:\n";
    outfile << indent << std::right
            << std::setw(label_width);
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width, '-') << '\n';
    outfile << indent << std::setw(label_width)
            << std::scientific << std::setprecision(4);
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << solver.sresidual1[i];
    }
    outfile << "\n\n";
}

void Output::write_inner_iteration(std::ofstream& outfile,
        coreData& core,    
        solverData& solver){

    const std::string indent = "    ";
    constexpr int cell_width = 14;
    const int table_width = cell_width * 2;
    const int niter = solver.ilinfphi.size();

    outfile << " Fixed source iteration summary:\n";
    outfile << "    Total number of iterations: " << niter << '\n';
    outfile << indent << std::left
            << std::setw(9) << "Iteration" << std::right
            << std::setw(cell_width) << "rho flux"
            << std::setw(cell_width) << "flux error"
            << '\n';
    outfile << indent << std::string(9+table_width, '-') << '\n';
    for (int i = 0; i < niter ; ++i) {
            outfile << indent << std::setw(9) << i+1
                << std::setw(cell_width) << std::scientific << std::setprecision(6) << solver.ispecrad[i]
                << std::setw(cell_width) << std::scientific << std::setprecision(6) << solver.ilinfphi[i]
                << '\n';
    }
    double stopcrit1 = core.epsflx * (1.0 / solver.ispecrad.back() - 1.0);
    outfile << "   Final flux criterion: " << std::scientific << std::setprecision(5)
                          << std::setw(16) << stopcrit1 << "\n"; 
    outfile << "\n\n";

}

void Output::write_outer_iteration(std::ofstream& outfile,
        coreData& core,
        solverData& solver){

    const std::string indent = "    ";
    constexpr int cell_width = 14;
    const int table_width = cell_width * 5;
    const int niter = solver.olinfphi.size();
       
    outfile << " Power iteration summary:\n";
    outfile << "    Total number of iterations: " << niter << '\n';
    outfile << indent << std::left
            << std::setw(9) << "Iteration" << std::right
            << std::setw(cell_width) << "keff"
            << std::setw(cell_width) << "keff error"
            << std::setw(cell_width) << "flux error"
            << std::setw(cell_width) << "rho flux"
            << std::setw(cell_width) << "fsource error"
            << '\n';
    outfile << indent << std::string(9+table_width, '-') << '\n';
    for (int i = 0; i < niter ; ++i) {
            outfile <<  std::setw(9) << i+1
                << std::setw(cell_width) << std::scientific << std::setprecision(6) << solver.keff[i]
                << std::setw(cell_width) << std::scientific << std::setprecision(6) << solver.linfkeff[i]
                << std::setw(cell_width) << std::scientific << std::setprecision(6) << solver.olinfphi[i]
                << std::setw(cell_width) << std::scientific << std::setprecision(6) << solver.ospecrad[i]
                << std::setw(cell_width) << std::scientific << std::setprecision(6) << solver.linffsrc[i]
                << '\n';
    }
    double stopcrit1 = core.epsk * (1.0 / solver.specradk.back() - 1.0);
    double stopcrit2 = core.epsflx * (1.0 / solver.ospecrad.back() - 1.0);
    outfile << "   Final keff criterion: " << std::scientific << std::setprecision(5)
                          << std::setw(16) << stopcrit1 << "\n";
    outfile << "   Final flux criterion: " << std::scientific << std::setprecision(5)
                          << std::setw(16) << stopcrit2 << "\n";                          
    outfile << "\n\n";

}