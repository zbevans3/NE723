#include "outer_iteration.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>

#include "source_iteration.h"

void OuterIteration::solve(coreData &core, quadratureData &quadrature, std::vector<materialData> &materials,
                           meshData &mesh, solverData &solver, QDData &QD, QDMGData &QDMG, QD1GData &QD1G) {

    double rhok = 1.0;
    double rhoflx = 1.0;
    double kold = 0.0;
    double keps = 10.0;
    double kepslm2 = 11.0;
    double maxepsflx = 10.0;
    double maxepsflxlm2 = 11.0;
    double sfistotnew = 1.0;
    double sfistotold = 0.0;
    double stopkcrit = 1.0e6;
    double stopflx = 1.0e6;
    double maxsfeps = 1.0e6;
    const int ngroups = core.ngroups;
    std::vector<std::vector<double>> phiold; // previous group scalar flux [g][i]
    std::vector<double> sfisnew;             // fission rate in cell i, summed over groups
    std::vector<double> sfisold;
    std::vector<std::vector<double>> epsflx; // group scalar flux change [g][i]
    std::vector<double> sfiseps;

    phiold.assign(ngroups, std::vector<double>(mesh.ncells, 1.0));
    epsflx.assign(ngroups, std::vector<double>(mesh.ncells, 0.0));
    sfisnew.assign(mesh.ncells, 1.0);
    sfisold.assign(mesh.ncells, 1.0);
    sfiseps.assign(mesh.ncells, 0.0);

    double invlambda = solver.lambda;
    SourceIteration inner;

    if (core.search != "critical_size")
        std::cout << "starting outer iteration\n";
    solver.k_outer = 0;

    if (core.otdebug) {
        std::cout << std::right << std::setw(20) << "Outer iteration" << std::setw(16) << "keff" << std::setw(16)
                  << "keff error" << std::setw(16) << "flux error" << std::setw(16) << "fsource error" << '\n';
    }

    int lout = 0;
    while (lout < core.maxout) {
        lout++;

        // update l-2 values
        kepslm2 = keps;
        maxepsflxlm2 = maxepsflx;

        kold = solver.lambda;
        phiold = solver.phibar;
        sfistotold = sfistotnew;
        sfisold = sfisnew;

        invlambda = 1.0 / solver.lambda;

        // do inner iteration
        inner.solve(core, quadrature, materials, mesh, solver, QD, QDMG, QD1G, invlambda);

        // compute new keff from new fission source
        sfistotnew = 0.0;
        maxsfeps = 0.0;
        maxepsflx = 0.0;
        for (int i = 0; i < mesh.ncells; i++) {
            // fission rate in cell i, summed over groups
            sfisnew[i] = 0.0;
            for (int g = 0; g < ngroups; g++) {
                sfisnew[i] +=
                    materials[mesh.mapmat[i]].nu[g] * materials[mesh.mapmat[i]].xsfiss[g] * solver.phibar[g][i];
                epsflx[g][i] = std::fabs(solver.phibar[g][i] - phiold[g][i]);
                maxepsflx = std::max(maxepsflx, epsflx[g][i]);
            }
            sfistotnew += sfisnew[i] * mesh.hx[i];
            sfiseps[i] = std::fabs(sfisnew[i] - sfisold[i]) / std::max(std::fabs(sfisold[i]), 1.0e-14);
            maxsfeps = std::max(maxsfeps, sfiseps[i]);
        }
        rhoflx = maxepsflx / maxepsflxlm2;

        solver.lambda = kold * (sfistotnew / sfistotold);
        solver.keff.push_back(solver.lambda);

        keps = fabs(solver.lambda - kold);
        rhok = keps / kepslm2;

        solver.linfkeff.push_back(keps);
        solver.specradk.push_back(rhok);
        solver.olinfphi.push_back(maxepsflx);
        solver.ospecrad.push_back(rhoflx);
        solver.linffsrc.push_back(maxsfeps);

        if (core.otdebug) {
            std::cout << std::right << std::setw(20) << lout << std::setw(16) << std::fixed << std::setprecision(8)
                      << solver.lambda << std::setw(16) << std::scientific << std::setprecision(6) << keps
                      << std::setw(16) << std::scientific << std::setprecision(6) << maxepsflx << std::setw(16)
                      << std::scientific << std::setprecision(6) << maxsfeps << '\n';
        }

        // stop if the solution blew up (inf or nan) instead of treating it as converged
        if (!std::isfinite(solver.lambda) || !std::isfinite(maxepsflx)) {
            throw std::runtime_error("Outer iteration diverged (keff or flux is inf/nan) at iteration " +
                                     std::to_string(lout));
        }

        // compute stopping criterion
        stopkcrit = core.epsk * (1.0 / rhok - 1.0);
        stopflx = core.epsflx * (1.0 / rhoflx - 1.0);
        if (keps < stopkcrit && maxepsflx < stopflx && maxsfeps < core.epspow) {
            if (core.search != "critical_size") {
                std::cout << "Outer iteration convergence reached in " << lout << " iterations.\n";
                std::cout << "Final convergence criterion: \n";
                std::cout << "  keff" << std::setw(16) << std::scientific << std::setprecision(6) << stopkcrit << "\n";
                std::cout << "  flux" << std::setw(16) << std::scientific << std::setprecision(6) << stopflx << "\n";
                std::cout << "\n";
            }
            break;
        }
    }

    // normalize so that (1/H) * sum_g sum_i phi_g,i * dx_i = 1
    const double H = mesh.xloc.back() - mesh.xloc.front();
    double phitot = 0.0;
    for (int g = 0; g < ngroups; g++) {
        for (int i = 0; i < mesh.ncells; i++) {
            phitot += solver.phibar[g][i] * mesh.hx[i];
        }
    }
    const double norm = H / phitot;

    // scale every group flux and current by the same factor
    for (int g = 0; g < ngroups; g++) {
        for (int i = 0; i < mesh.ncells; i++) {
            solver.phibar[g][i] *= norm;
            solver.phihat[g][i] *= norm;
            solver.jbar[g][i] *= norm;
            solver.jhat[g][i] *= norm;
        }
        for (int i = 0; i < mesh.nedges; i++) {
            solver.phie[g][i] *= norm;
            solver.je[g][i] *= norm;
        }
    }

    solver.k_outer += lout;
}