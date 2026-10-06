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
        outfile << "    material " << std::setw(6) << materials[i].id << std::setw(8) 
                << std::fixed << std::setprecision(4) << materials[i].xstot << std::setw(8)
                << std::fixed << std::setprecision(4) << materials[i].xsscat << std::setw(8)
                << std::fixed << std::setprecision(4) << materials[i].xsfiss << std::setw(8)
                << std::fixed << std::setprecision(4) << materials[i].nu << std::setw(8)
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
    outfile << "\n\n\n";
}

void Output::write_mesh(
    std::ofstream& outfile,
    coreData& core,
    meshData& mesh)
{
    double center = 0.0;

    outfile << " Mesh information:\n";
    outfile << "    Number of cells: " << mesh.ncells << '\n';
    outfile << "    Number of edges: " << mesh.nedges << "\n\n";

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
        center = 0.5 * (mesh.xloc[i] + mesh.xloc[i + 1]);
        outfile << std::setw(8)  << i + 1
                << std::setw(10) << mesh.mapmat[i]
                << std::setw(15) << mesh.xloc[i]
                << std::setw(15) << mesh.xloc[i + 1]
                << std::setw(15) << center
                << std::setw(15) << mesh.hx[i]
                << '\n';
    }

    outfile << "\n Edge information:\n";
    outfile << std::setw(8)  << "Edge" << std::setw(15) << "Location" << std::setw(10) << "ibc" << '\n';
    outfile << "    " << std::string(29, '-') << '\n';
    for (int i = 0; i < mesh.nedges; ++i) {
        outfile << std::setw(8)  << i + 1 << std::setw(15) << mesh.xloc[i] << std::setw(10) << mesh.mapbc[i] << '\n';
    }
    outfile << '\n';
}