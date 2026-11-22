#include "write_output.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <sstream>

void Output::write_header(std::ofstream &outfile) {
    outfile << "===================================================================\n";
    outfile << "=                            NE723 Project 3                      =\n";
    outfile << "=                            Zachary Bevans                       =\n";
    outfile << "=                              12/01/2026                         =\n";
    outfile << "===================================================================\n";
    outfile << "Code to solve the steady state multigroup neutron transport equation\n";
    outfile << " in 1D slab geometry using the discrete ordinates method \n";
    outfile << " and Quasi-Diffusion (QD) acceleration.\n";
    outfile << "===================================================================\n\n";
}

void Output::echo_input(std::ofstream &outfile, coreData &core, quadratureData &quadrature,
                        std::vector<materialData> &materials) {

    outfile << " Input values:\n";
    outfile << "  [CORE]\n";
    outfile << std::left;
    outfile << "    " << std::setw(16) << "search" << core.search << "\n";
    outfile << "    " << std::setw(16) << "acceleration" << core.acceleration << "\n";
    outfile << "    " << std::setw(16) << "ngroups" << core.ngroups << "\n";
    outfile << "    " << std::setw(16) << "itdebug" << core.itdebug << "\n";
    outfile << "    " << std::setw(16) << "epsk" << core.epsk << "\n";
    outfile << "    " << std::setw(16) << "epsflx" << core.epsflx << "\n";
    outfile << "    " << std::setw(16) << "maxin" << core.maxin << "\n";
    outfile << "    " << std::setw(16) << "maxout" << core.maxout << "\n";
    outfile << "    " << std::setw(16) << "nregions" << core.nregions << "\n";
    outfile << "    " << std::setw(16) << "nmat" << core.nmat << "\n";
    outfile << "    " << std::setw(16) << "bcleft" << std::setw(16) << core.bcleft << std::setw(6) << std::fixed
            << std::setprecision(3) << core.psileft << "\n";
    outfile << "    " << std::setw(16) << "bcright" << std::setw(16) << core.bcright << std::setw(5) << std::fixed
            << std::setprecision(3) << core.psiright << "\n";
    outfile << std::right;
    outfile << "    nx ";
    for (int i = 0; i < core.nregions; i++) {
        outfile << std::setw(6) << core.nx[i];
    }
    outfile << "\n";
    outfile << "    xedge ";
    for (int i = 0; i < core.nregions; i++) {
        outfile << std::setw(8) << std::fixed << std::setprecision(3) << core.xedge[i];
    }
    outfile << "\n\n";

    outfile << "  [MATERIAL]\n";
    // label column is 14 wide, each group column is 13 wide
    int _ngroups = core.ngroups;
    outfile << std::scientific << std::setprecision(7);
    for (int i = 0; i < core.nmat; i++) {
        outfile << "    material " << materials[i].id << "\n";
        outfile << "    " << std::left << std::setw(14) << "group" << std::right;
        for (int g = 0; g < _ngroups; g++) {
            outfile << std::setw(16) << g + 1;
        }
        outfile << "\n";
        outfile << "    " << std::left << std::setw(14) << "total" << std::right;
        for (int g = 0; g < _ngroups; g++) {
            outfile << std::setw(16) << materials[i].xstot[g];
        }
        outfile << "\n";
        outfile << "    " << std::left << std::setw(14) << "absorption" << std::right;
        for (int g = 0; g < _ngroups; g++) {
            outfile << std::setw(16) << materials[i].xsabs[g];
        }
        outfile << "\n";
        outfile << "    " << std::left << std::setw(14) << "fission" << std::right;
        for (int g = 0; g < _ngroups; g++) {
            outfile << std::setw(16) << materials[i].xsfiss[g];
        }
        outfile << "\n";
        outfile << "    " << std::left << std::setw(14) << "nu" << std::right;
        for (int g = 0; g < _ngroups; g++) {
            outfile << std::setw(16) << materials[i].nu[g];
        }
        outfile << "\n";
        outfile << "    " << std::left << std::setw(14) << "chi" << std::right;
        for (int g = 0; g < _ngroups; g++) {
            outfile << std::setw(16) << materials[i].xschi[g];
        }
        outfile << "\n";
        outfile << "    " << std::left << std::setw(14) << "fixed_source" << std::right;
        for (int g = 0; g < _ngroups; g++) {
            outfile << std::setw(16) << materials[i].q[g];
        }
        outfile << "\n";
        outfile << "    scattering matrix (row = from group, column = to group)\n";
        outfile << "    " << std::left << std::setw(14) << "from \\ to" << std::right;
        for (int g = 0; g < _ngroups; g++) {
            outfile << std::setw(16) << g + 1;
        }
        outfile << "\n";
        for (int g = 0; g < _ngroups; g++) {
            outfile << "    " << std::left << std::setw(14) << g + 1 << std::right;
            for (int gp = 0; gp < _ngroups; gp++) {
                outfile << std::setw(16) << materials[i].xsscatmat[g][gp];
            }
            outfile << "\n";
        }
        outfile << "\n";
    }

    outfile << "  [QUADRATURE]\n";
    outfile << "    order " << std::setw(6) << std::right << quadrature.order << "\n";
    outfile << "    angle ";
    for (int i = 0; i < quadrature.order; i++) {
        outfile << std::setw(16) << std::right << std::scientific << std::setprecision(7) << quadrature.angle[i];
    }
    outfile << "\n";
    outfile << "    weight ";
    for (int i = 0; i < quadrature.order; i++) {
        outfile << std::setw(16) << std::right << std::scientific << std::setprecision(7) << quadrature.weight[i];
    }
    outfile << "\n";
    outfile << "    sum of weights " << std::setw(16) << std::right << std::scientific << std::setprecision(7)
            << std::accumulate(quadrature.weight.begin(), quadrature.weight.end(), 0.0) << "\n";
    double wsum = 0.0;
    for (int i = 0; i < quadrature.order; i++) {
        wsum += quadrature.angle[i] * quadrature.weight[i];
    }
    outfile << "    sum of weight*angle " << std::setw(16) << std::right << std::scientific << std::setprecision(7)
            << wsum << "\n";
    outfile << "\n\n";
}

void Output::write_mesh(std::ofstream &outfile, meshData &mesh, std::vector<materialData> &materials) {

    outfile << " Cell information:\n";
    outfile << "    Number of cells: " << mesh.ncells << '\n';

    outfile << std::right << std::setw(8) << "Cell" << std::setw(10) << "Material" << std::setw(16) << "Left edge"
            << std::setw(16) << "Right edge" << std::setw(16) << "Center" << std::setw(16) << "Width" << '\n';

    outfile << "    " << std::string(82, '-') << '\n';
    outfile << std::scientific << std::setprecision(7);

    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(8) << i + 1 << std::setw(10) << materials[mesh.mapmat[i]].id << std::setw(16)
                << mesh.xloc[i] << std::setw(16) << mesh.xloc[i + 1] << std::setw(16) << mesh.xcenter[i]
                << std::setw(16) << mesh.hx[i] << '\n';
    }

    outfile << "\n Edge information:\n";
    outfile << "    Number of edges: " << mesh.nedges << "\n";
    outfile << std::setw(8) << "Edge" << std::setw(16) << "Location" << std::setw(10) << "ibc" << '\n';
    outfile << "    " << std::string(34, '-') << '\n';
    for (int i = 0; i < mesh.nedges; ++i) {
        outfile << std::setw(8) << i + 1 << std::setw(16) << mesh.xloc[i] << std::setw(10) << mesh.mapbc[i] << '\n';
    }
    outfile << '\n';
}

void Output::write_angular_flux(std::ofstream &outfile, meshData &mesh, quadratureData &quadrature,
                                solverData &solver) {
    const std::string indent = "    ";
    constexpr int label_width = 16;
    constexpr int cell_width = 16;
    const int table_width = label_width + cell_width * mesh.ncells;
    const int table_width_e = label_width + cell_width * mesh.nedges;
    const int ngroups = solver.psibar.size();

    for (int g = 0; g < ngroups; ++g) {
        outfile << " Cell average angular flux solution, group " << g + 1 << ":\n";
        outfile << indent << std::right << std::setw(label_width) << "Angle";
        for (int i = 0; i < mesh.ncells; ++i) {
            outfile << std::setw(cell_width) << i + 1;
        }
        outfile << '\n';
        outfile << indent << std::string(table_width, '-') << '\n';
        outfile << std::scientific << std::setprecision(7);
        for (int u = 0; u < quadrature.order; ++u) {
            outfile << indent << std::scientific << std::setprecision(7) << std::setw(label_width)
                    << quadrature.angle[u] << std::scientific << std::setprecision(7);
            for (int i = 0; i < mesh.ncells; ++i) {
                outfile << std::setw(cell_width) << solver.psibar[g][u][i];
            }
            outfile << '\n';
        }
        outfile << '\n';
    }

    for (int g = 0; g < ngroups; ++g) {
        outfile << " First Legendre spatial moment of the angular flux, group " << g + 1 << ":\n";
        outfile << indent << std::right << std::setw(label_width) << "Angle";
        for (int i = 0; i < mesh.ncells; ++i) {
            outfile << std::setw(cell_width) << i + 1;
        }
        outfile << '\n';
        outfile << indent << std::string(table_width, '-') << '\n';
        outfile << std::scientific << std::setprecision(7);
        for (int u = 0; u < quadrature.order; ++u) {
            outfile << indent << std::scientific << std::setprecision(7) << std::setw(label_width)
                    << quadrature.angle[u] << std::scientific << std::setprecision(7);
            for (int i = 0; i < mesh.ncells; ++i) {
                outfile << std::setw(cell_width) << solver.psihat[g][u][i];
            }
            outfile << '\n';
        }
        outfile << '\n';
    }

    for (int g = 0; g < ngroups; ++g) {
        outfile << " Cell edge angular flux solution, group " << g + 1 << ":\n";
        outfile << indent << std::right << std::setw(label_width) << "Angle";
        for (int i = 0; i < mesh.nedges; ++i) {
            outfile << std::setw(cell_width) << i + 1;
        }
        outfile << '\n';
        outfile << indent << std::string(table_width_e, '-') << '\n';
        outfile << std::scientific << std::setprecision(7);
        for (int u = 0; u < quadrature.order; ++u) {
            outfile << indent << std::scientific << std::setprecision(7) << std::setw(label_width)
                    << quadrature.angle[u] << std::scientific << std::setprecision(7);
            for (int i = 0; i < mesh.nedges; ++i) {
                outfile << std::setw(cell_width) << solver.psie[g][u][i];
            }
            outfile << '\n';
        }
        outfile << '\n';
    }
}

void Output::write_scalar_flux(std::ofstream &outfile, meshData &mesh, solverData &solver) {
    const std::string indent = "    ";
    constexpr int label_width = 8;
    constexpr int cell_width = 16;
    const int table_width = label_width + cell_width * mesh.ncells;
    const int table_width_e = label_width + cell_width * mesh.nedges;
    const int ngroups = solver.phibar.size();

    outfile << " Cell average scalar flux solution:\n";
    outfile << indent << std::right << std::setw(label_width) << "Group";
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width, '-') << '\n';
    outfile << std::scientific << std::setprecision(7);
    for (int g = 0; g < ngroups; ++g) {
        outfile << indent << std::setw(label_width) << g + 1;
        for (int i = 0; i < mesh.ncells; ++i) {
            outfile << std::setw(cell_width) << solver.phibar[g][i];
        }
        outfile << '\n';
    }
    outfile << '\n';

    outfile << " First Legendre spatial moment of the scalar flux:\n";
    outfile << indent << std::right << std::setw(label_width) << "Group";
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width, '-') << '\n';
    outfile << std::scientific << std::setprecision(7);
    for (int g = 0; g < ngroups; ++g) {
        outfile << indent << std::setw(label_width) << g + 1;
        for (int i = 0; i < mesh.ncells; ++i) {
            outfile << std::setw(cell_width) << solver.phihat[g][i];
        }
        outfile << '\n';
    }
    outfile << '\n';

    outfile << " Cell edge scalar flux solution:\n";
    outfile << indent << std::right << std::setw(label_width) << "Group";
    for (int i = 0; i < mesh.nedges; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width_e, '-') << '\n';
    outfile << std::scientific << std::setprecision(7);
    for (int g = 0; g < ngroups; ++g) {
        outfile << indent << std::setw(label_width) << g + 1;
        for (int i = 0; i < mesh.nedges; ++i) {
            outfile << std::setw(cell_width) << solver.phie[g][i];
        }
        outfile << '\n';
    }
    outfile << '\n';
}

void Output::write_current(std::ofstream &outfile, meshData &mesh, solverData &solver) {
    const std::string indent = "    ";
    constexpr int label_width = 8;
    constexpr int cell_width = 16;
    const int table_width = label_width + cell_width * mesh.ncells;
    const int table_width_e = label_width + cell_width * mesh.nedges;
    const int ngroups = solver.jbar.size();

    outfile << " Cell average current solution:\n";
    outfile << indent << std::right << std::setw(label_width) << "Group";
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width, '-') << '\n';
    outfile << std::scientific << std::setprecision(7);
    for (int g = 0; g < ngroups; ++g) {
        outfile << indent << std::setw(label_width) << g + 1;
        for (int i = 0; i < mesh.ncells; ++i) {
            outfile << std::setw(cell_width) << solver.jbar[g][i];
        }
        outfile << '\n';
    }
    outfile << '\n';

    outfile << " First Legendre spatial moment of the current:\n";
    outfile << indent << std::right << std::setw(label_width) << "Group";
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width, '-') << '\n';
    outfile << std::scientific << std::setprecision(7);
    for (int g = 0; g < ngroups; ++g) {
        outfile << indent << std::setw(label_width) << g + 1;
        for (int i = 0; i < mesh.ncells; ++i) {
            outfile << std::setw(cell_width) << solver.jhat[g][i];
        }
        outfile << '\n';
    }
    outfile << '\n';

    outfile << " Cell edge current solution:\n";
    outfile << indent << std::right << std::setw(label_width) << "Group";
    for (int i = 0; i < mesh.nedges; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width_e, '-') << '\n';
    outfile << std::scientific << std::setprecision(7);
    for (int g = 0; g < ngroups; ++g) {
        outfile << indent << std::setw(label_width) << g + 1;
        for (int i = 0; i < mesh.nedges; ++i) {
            outfile << std::setw(cell_width) << solver.je[g][i];
        }
        outfile << '\n';
    }
    outfile << '\n';
}

void Output::write_scalar_flux_LO(std::ofstream &outfile, meshData &mesh, QDData &QD) {
    const std::string indent = "    ";
    constexpr int label_width = 8;
    constexpr int cell_width = 16;
    const int table_width = cell_width * mesh.ncells;
    const int table_width_e = cell_width * mesh.nedges;

    outfile << " Cell average LO scalar flux solution:\n";
    outfile << indent << std::right << std::setw(label_width);
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width, '-') << '\n';
    outfile << indent << std::setw(label_width) << std::scientific << std::setprecision(7);
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << QD.phibarLO[i];
    }
    outfile << "\n\n";

    outfile << " Cell edge LO scalar flux solution:\n";
    outfile << indent << std::right << std::setw(label_width);
    for (int i = 0; i < mesh.nedges; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width_e, '-') << '\n';
    outfile << indent << std::setw(label_width) << std::scientific << std::setprecision(7);
    for (int i = 0; i < mesh.nedges; ++i) {
        outfile << std::setw(cell_width) << QD.phieLO[i];
    }
    outfile << "\n\n";
}

void Output::write_current_LO(std::ofstream &outfile, meshData &mesh, QDData &QD) {
    const std::string indent = "    ";
    constexpr int label_width = 8;
    constexpr int cell_width = 16;
    const int table_width_e = cell_width * mesh.nedges;

    outfile << " Cell edge LO current solution:\n";
    outfile << indent << std::right << std::setw(label_width);
    for (int i = 0; i < mesh.nedges; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width_e, '-') << '\n';
    outfile << indent << std::setw(label_width) << std::scientific << std::setprecision(7);
    for (int i = 0; i < mesh.nedges; ++i) {
        outfile << std::setw(cell_width) << QD.jeLO[i];
    }
    outfile << "\n\n";
}

void Output::write_QD_data(std::ofstream &outfile, meshData &mesh, QDData &QD) {

    const std::string indent = "    ";
    constexpr int label_width = 8;
    constexpr int cell_width = 16;
    const int table_width = cell_width * mesh.ncells;
    const int table_width_e = cell_width * mesh.nedges;

    outfile << " QD data:\n";
    outfile << std::left;
    outfile << std::scientific << std::setprecision(7);
    outfile << "  " << std::setw(24) << "Left boundary factor" << QD.CbfL << "\n";
    outfile << "  " << std::setw(24) << "Right boundary factor" << QD.CbfR << "\n";
    outfile << "\n";

    outfile << " Cell average Eddington factor:\n";
    outfile << indent << std::right << std::setw(label_width);
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width, '-') << '\n';
    outfile << indent << std::setw(label_width) << std::scientific << std::setprecision(7);
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << QD.Eddbar[i];
    }
    outfile << "\n\n";

    outfile << " Cell edge Eddington factor:\n";
    outfile << indent << std::right << std::setw(label_width);
    for (int i = 0; i < mesh.nedges; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width_e, '-') << '\n';
    outfile << indent << std::setw(label_width) << std::scientific << std::setprecision(7);
    for (int i = 0; i < mesh.nedges; ++i) {
        outfile << std::setw(cell_width) << QD.Edde[i];
    }
    outfile << "\n\n";

    outfile << " Metric of transport effects:\n";
    outfile << indent << std::right << std::setw(label_width);
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width, '-') << '\n';
    outfile << indent << std::setw(label_width) << std::scientific << std::setprecision(7);
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << QD.mLO[i];
    }
    outfile << "\n\n";

    outfile << " Derivative based metric of transport effects:\n";
    outfile << indent << std::right << std::setw(label_width);
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width, '-') << '\n';
    outfile << indent << std::setw(label_width) << std::scientific << std::setprecision(7);
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << QD.mprimeLO[i];
    }
    outfile << "\n\n";
}

void Output::write_QDMG_data(std::ofstream &outfile, meshData &mesh, QDMGData &QDMG) {

    const std::string indent = "    ";
    constexpr int label_width = 8;
    constexpr int cell_width = 16;
    const int table_width = label_width + cell_width * mesh.ncells;
    const int table_width_e = label_width + cell_width * mesh.nedges;
    const int ngroups = QDMG.phibar.size();

    outfile << " Multigroup QD data:\n";
    outfile << indent << std::right << std::setw(label_width) << "Group" << std::setw(cell_width) << "Left factor"
            << std::setw(cell_width) << "Right factor" << '\n';
    outfile << indent << std::string(label_width + 2 * cell_width, '-') << '\n';
    outfile << std::scientific << std::setprecision(7);
    for (int g = 0; g < ngroups; ++g) {
        outfile << indent << std::setw(label_width) << g + 1 << std::setw(cell_width) << QDMG.CbfL[g]
                << std::setw(cell_width) << QDMG.CbfR[g] << '\n';
    }
    outfile << '\n';

    outfile << " Cell average group Eddington factor:\n";
    outfile << indent << std::right << std::setw(label_width) << "Group";
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width, '-') << '\n';
    outfile << std::scientific << std::setprecision(7);
    for (int g = 0; g < ngroups; ++g) {
        outfile << indent << std::setw(label_width) << g + 1;
        for (int i = 0; i < mesh.ncells; ++i) {
            outfile << std::setw(cell_width) << QDMG.Eddbar[g][i];
        }
        outfile << '\n';
    }
    outfile << '\n';

    outfile << " Cell edge group Eddington factor:\n";
    outfile << indent << std::right << std::setw(label_width) << "Group";
    for (int i = 0; i < mesh.nedges; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width_e, '-') << '\n';
    outfile << std::scientific << std::setprecision(7);
    for (int g = 0; g < ngroups; ++g) {
        outfile << indent << std::setw(label_width) << g + 1;
        for (int i = 0; i < mesh.nedges; ++i) {
            outfile << std::setw(cell_width) << QDMG.Edde[g][i];
        }
        outfile << '\n';
    }
    outfile << '\n';

    outfile << " Cell average group LO scalar flux solution:\n";
    outfile << indent << std::right << std::setw(label_width) << "Group";
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width, '-') << '\n';
    outfile << std::scientific << std::setprecision(7);
    for (int g = 0; g < ngroups; ++g) {
        outfile << indent << std::setw(label_width) << g + 1;
        for (int i = 0; i < mesh.ncells; ++i) {
            outfile << std::setw(cell_width) << QDMG.phibar[g][i];
        }
        outfile << '\n';
    }
    outfile << '\n';

    outfile << " Cell edge group LO scalar flux solution:\n";
    outfile << indent << std::right << std::setw(label_width) << "Group";
    for (int i = 0; i < mesh.nedges; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width_e, '-') << '\n';
    outfile << std::scientific << std::setprecision(7);
    for (int g = 0; g < ngroups; ++g) {
        outfile << indent << std::setw(label_width) << g + 1;
        for (int i = 0; i < mesh.nedges; ++i) {
            outfile << std::setw(cell_width) << QDMG.phie[g][i];
        }
        outfile << '\n';
    }
    outfile << '\n';

    outfile << " Cell edge group LO current solution:\n";
    outfile << indent << std::right << std::setw(label_width) << "Group";
    for (int i = 0; i < mesh.nedges; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width_e, '-') << '\n';
    outfile << std::scientific << std::setprecision(7);
    for (int g = 0; g < ngroups; ++g) {
        outfile << indent << std::setw(label_width) << g + 1;
        for (int i = 0; i < mesh.nedges; ++i) {
            outfile << std::setw(cell_width) << QDMG.je[g][i];
        }
        outfile << '\n';
    }
    outfile << '\n';
}

void Output::write_QD1G_data(std::ofstream &outfile, meshData &mesh, QD1GData &QD1G) {

    const std::string indent = "    ";
    constexpr int label_width = 8;
    constexpr int cell_width = 16;
    const int table_width = cell_width * mesh.ncells;
    const int table_width_e = cell_width * mesh.nedges;

    outfile << " One-group QD data:\n";
    outfile << std::left;
    outfile << std::scientific << std::setprecision(7);
    outfile << "  " << std::setw(24) << "One-group eigenvalue" << QD1G.k << "\n";
    outfile << "  " << std::setw(24) << "Left boundary factor" << QD1G.CbfL << "\n";
    outfile << "  " << std::setw(24) << "Right boundary factor" << QD1G.CbfR << "\n";
    outfile << "\n";

    outfile << " Cell average one-group LO scalar flux solution:\n";
    outfile << indent << std::right << std::setw(label_width);
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width, '-') << '\n';
    outfile << indent << std::setw(label_width) << std::scientific << std::setprecision(7);
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << QD1G.phibar[i];
    }
    outfile << "\n\n";

    outfile << " Cell edge one-group LO scalar flux solution:\n";
    outfile << indent << std::right << std::setw(label_width);
    for (int i = 0; i < mesh.nedges; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width_e, '-') << '\n';
    outfile << indent << std::setw(label_width) << std::scientific << std::setprecision(7);
    for (int i = 0; i < mesh.nedges; ++i) {
        outfile << std::setw(cell_width) << QD1G.phie[i];
    }
    outfile << "\n\n";

    outfile << " Cell edge one-group LO current solution:\n";
    outfile << indent << std::right << std::setw(label_width);
    for (int i = 0; i < mesh.nedges; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width_e, '-') << '\n';
    outfile << indent << std::setw(label_width) << std::scientific << std::setprecision(7);
    for (int i = 0; i < mesh.nedges; ++i) {
        outfile << std::setw(cell_width) << QD1G.je[i];
    }
    outfile << "\n\n";
}

void Output::write_residuals(std::ofstream &outfile, meshData &mesh, solverData &solver) {
    const std::string indent = "    ";
    constexpr int label_width = 8;
    constexpr int cell_width = 16;
    const int table_width = label_width + cell_width * mesh.ncells;
    const int ngroups = solver.sresidual.size();

    outfile << " Cell average scheme residual:\n";
    outfile << indent << std::right << std::setw(label_width) << "Group";
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width, '-') << '\n';
    outfile << std::scientific << std::setprecision(7);
    for (int g = 0; g < ngroups; ++g) {
        outfile << indent << std::setw(label_width) << g + 1;
        for (int i = 0; i < mesh.ncells; ++i) {
            outfile << std::setw(cell_width) << solver.sresidual[g][i];
        }
        outfile << '\n';
    }
    outfile << '\n';

    outfile << " Scheme residual for first Legendre spatial moment:\n";
    outfile << indent << std::right << std::setw(label_width) << "Group";
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width, '-') << '\n';
    outfile << std::scientific << std::setprecision(7);
    for (int g = 0; g < ngroups; ++g) {
        outfile << indent << std::setw(label_width) << g + 1;
        for (int i = 0; i < mesh.ncells; ++i) {
            outfile << std::setw(cell_width) << solver.sresidual1[g][i];
        }
        outfile << '\n';
    }
    outfile << '\n';

    outfile << " Cell average iterative residual:\n";
    outfile << indent << std::right << std::setw(label_width) << "Group";
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width, '-') << '\n';
    outfile << std::scientific << std::setprecision(7);
    for (int g = 0; g < ngroups; ++g) {
        outfile << indent << std::setw(label_width) << g + 1;
        for (int i = 0; i < mesh.ncells; ++i) {
            outfile << std::setw(cell_width) << solver.iresidual[g][i];
        }
        outfile << '\n';
    }
    outfile << '\n';

    outfile << " Iterative residual for first Legendre spatial moment:\n";
    outfile << indent << std::right << std::setw(label_width) << "Group";
    for (int i = 0; i < mesh.ncells; ++i) {
        outfile << std::setw(cell_width) << i + 1;
    }
    outfile << '\n';
    outfile << indent << std::string(table_width, '-') << '\n';
    outfile << std::scientific << std::setprecision(7);
    for (int g = 0; g < ngroups; ++g) {
        outfile << indent << std::setw(label_width) << g + 1;
        for (int i = 0; i < mesh.ncells; ++i) {
            outfile << std::setw(cell_width) << solver.iresidual1[g][i];
        }
        outfile << '\n';
    }
    outfile << '\n';
}

void Output::write_inner_iteration(std::ofstream &outfile, coreData &core, solverData &solver) {

    const std::string indent = "    ";
    constexpr int cell_width = 16;
    const int table_width = cell_width * 2;
    const int niter = solver.ilinfphi.size();

    outfile << " Fixed source iteration summary:\n";
    outfile << "    Total number of iterations: " << niter << '\n';
    outfile << indent << std::left << std::setw(9) << "Iteration" << std::right << std::setw(cell_width) << "rho flux"
            << std::setw(cell_width) << "flux error" << '\n';
    outfile << indent << std::string(9 + table_width, '-') << '\n';
    for (int i = 0; i < niter; ++i) {
        outfile << indent << std::setw(9) << i + 1 << std::setw(cell_width) << std::scientific << std::setprecision(7)
                << solver.ispecrad[i] << std::setw(cell_width) << std::scientific << std::setprecision(7)
                << solver.ilinfphi[i] << '\n';
    }
    double stopcrit1 = core.epsflx * (1.0 / solver.ispecrad.back() - 1.0);
    outfile << "   Final flux criterion: " << std::scientific << std::setprecision(7) << std::setw(16) << stopcrit1
            << "\n";
    outfile << "\n\n";
}

void Output::write_outer_iteration(std::ofstream &outfile, coreData &core, solverData &solver) {

    const std::string indent = "    ";
    constexpr int cell_width = 16;
    const int table_width = cell_width * 5;
    const int niter = solver.olinfphi.size();

    outfile << " Power iteration summary:\n";
    outfile << "    Total number of iterations: " << niter << '\n';
    outfile << "    Total number of inner iterations: " << solver.k_inner << '\n';
    outfile << indent << std::left << std::setw(9) << "Iteration" << std::right << std::setw(cell_width) << "keff"
            << std::setw(cell_width) << "keff error" << std::setw(cell_width) << "flux error" << std::setw(cell_width)
            << "rho flux" << std::setw(cell_width) << "fsource error" << '\n';
    outfile << indent << std::string(9 + table_width, '-') << '\n';
    for (int i = 0; i < niter; ++i) {
        outfile << std::setw(9) << i + 1 << std::setw(cell_width) << std::scientific << std::setprecision(7)
                << solver.keff[i] << std::setw(cell_width) << std::scientific << std::setprecision(7)
                << solver.linfkeff[i] << std::setw(cell_width) << std::scientific << std::setprecision(7)
                << solver.olinfphi[i] << std::setw(cell_width) << std::scientific << std::setprecision(7)
                << solver.ospecrad[i] << std::setw(cell_width) << std::scientific << std::setprecision(7)
                << solver.linffsrc[i] << '\n';
    }
    double stopcrit1 = core.epsk * (1.0 / solver.specradk.back() - 1.0);
    double stopcrit2 = core.epsflx * (1.0 / solver.ospecrad.back() - 1.0);
    outfile << "   Final keff criterion: " << std::scientific << std::setprecision(7) << std::setw(16) << stopcrit1
            << "\n";
    outfile << "   Final flux criterion: " << std::scientific << std::setprecision(7) << std::setw(16) << stopcrit2
            << "\n";
    outfile << "\n\n";
}
