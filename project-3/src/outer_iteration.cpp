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
    double phimax = 1.0e6;
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

        if (core.acceleration == "MLQD"){
            // solve 1G QD system

            // first, compute flux-weighted one group cross sections and C_L/C_R
            QD1G.xsabs.assign(mesh.ncells, 0.0);
            QD1G.xsfiss.assign(mesh.ncells, 0.0);
            QD1G.nu.assign(mesh.ncells, 0.0);
            for (int i = 0 ; i < mesh.ncells ; i++){
                const materialData &mat = materials[mesh.mapmat[i]];
                double phisum = 0.0;
                double phifis = 0.0;         
                 for (int g = 0 ; g < ngroups ; g++){
                    QD1G.xsabs[i] += mat.xsabs[g] * QDMG.phibar[g][i];
                    QD1G.xsfiss[i] += mat.xsfiss[g] * QDMG.phibar[g][i];
                    QD1G.nu[i] += mat.nu[g] * mat.xsfiss[g] * QDMG.phibar[g][i];
                    phisum += QDMG.phibar[g][i];
                    phifis += mat.xsfiss[g] * QDMG.phibar[g][i];
                }
                QD1G.xsabs[i] /= phisum;
                QD1G.xsfiss[i] /= phisum;
                QD1G.nu[i] /= phifis;            
            }

            // Compute 1G boundary factors
            double phisumL = 0.0;
            double phisumR = 0.0;
            double csumL = 0.0;
            double csumR = 0.0;
            for (int g = 0 ; g < ngroups ; g++){
                csumL += QDMG.CbfL[g] * QDMG.phie[g].front();
                phisumL += QDMG.phie[g].front();
                csumR += QDMG.CbfR[g] * QDMG.phie[g].back();
                phisumR += QDMG.phie[g].back();
            }
            QD1G.CbfL = (mesh.mapbc.front() != 2) ? csumL / phisumL : 0.0;
            QD1G.CbfR = (mesh.mapbc.back() != 2) ? csumR / phisumR : 0.0;
           
            // compute Hm, Hp, Am, Ap
            QD1G.Hm.assign(mesh.nedges, 0.0);
            QD1G.Hp.assign(mesh.nedges, 0.0);
            QD1G.Am.assign(mesh.nedges, 0.0); 
            QD1G.Ap.assign(mesh.nedges, 0.0);
            for (int i = 0 ; i <= mesh.ncells ; i++){
                const int iL = (i == 0) ? 0 : i - 1;
                const int iR = (i == mesh.ncells) ? mesh.ncells - 1 : i;
                const materialData &matL = materials[mesh.mapmat[iL]];
                const materialData &matR = materials[mesh.mapmat[iR]];
                const double hL = (i == 0) ? 0.0 : mesh.hx[iL];
                const double hR = (i == mesh.ncells) ? 0.0 : mesh.hx[iR];
                double nsumm = 0.0;
                double dsumm = 0.0;
                double nsump = 0.0;
                double dsump = 0.0;
                double nsumam = 0.0;
                double nsumap = 0.0;
                for (int g = 0 ; g < ngroups ; g++){
                    const double sigmat = (matL.xstot[g] * hL + matR.xstot[g] * hR) / (hL + hR);
                    const double Edde = QDMG.Edde[g][i];
                    const double phiL = (i == 0) ? QDMG.phie[g][0] : QDMG.phibar[g][iL];
                    const double EddL = (i == 0) ? QDMG.Edde[g][0] : QDMG.Eddbar[g][iL];
                    const double phiR = (i == mesh.ncells) ? QDMG.phie[g][mesh.ncells] : QDMG.phibar[g][iR];
                    const double EddR = (i == mesh.ncells) ? QDMG.Edde[g][mesh.ncells] : QDMG.Eddbar[g][iR];
                    nsumm += EddL * phiL / sigmat;
                    dsumm += phiL;
                    nsump += EddR * phiR / sigmat;
                    dsump += phiR;
                    nsumam += matR.xstot[g] * hR * EddL * phiL / ((matL.xstot[g] * hL + matR.xstot[g] * hR) * Edde);
                    nsumap += matL.xstot[g] * hL * EddR * phiR / ((matL.xstot[g] * hL + matR.xstot[g] * hR) * Edde);
                }
                QD1G.Hm[i] = nsumm / dsumm;
                QD1G.Hp[i] = nsump / dsump;
                QD1G.Am[i] = nsumam / dsumm;
                QD1G.Ap[i] = nsumap / dsump;
            }

            // build one group system
            // we only solve for cell average values in linear system
            // then reconstruct edge fluxes/currents
            const int N = mesh.ncells;
            const int n = N;
            Al.assign(n, 0.0);
            Ad.assign(n, 0.0);
            Au.assign(n, 0.0);
            bvec.assign(n, 0.0);
            phivec.assign(n, 0.0);
            a1.assign(N, 0.0);
            a2.assign(N, 0.0);
            a3.assign(N, 0.0);
            b.assign(N, 0.0);

            // compute half-delta x
            hxe.assign(mesh.nedges, 0.0);
            for (int i = 0 ; i <= N ; i++){
                const double hL = (i == 0) ? 0.0 : mesh.hx[i - 1];
                const double hR = (i == N) ? 0.0 : mesh.hx[i];
                hxe[i] = 0.5 * (hL + hR);
            }
            
            // compute coupling coefficients
            for (int i = 0 ; i < N ; i++){
                a1[i] = -QD1G.Hm[i] / hxe[i];
                a2[i] = QD1G.Hm[i + 1] / hxe[i + 1] + QD1G.Hp[i] / hxe[i] + QD1G.xsabs[i] * mesh.hx[i];
                a3[i] = -QD1G.Hp[i + 1] / hxe[i + 1];
            }

            // left edge
            a1[0] = 0.0;
            a2[0] -= QD1G.Hp[0] / hxe[0];
            if (mesh.mapbc.front() != 2){
                a2[0] -= QD1G.CbfL * QD1G.Hp[0] / (QD1G.Hm[0] - QD1G.CbfL * hxe[0]);
            }

            // right edge
            a3[N - 1] = 0.0;
            a2[N - 1] -= QD1G.Hm[N] / hxe[N];
            if (mesh.mapbc.back() != 2){
                a2[N - 1] += QD1G.CbfR * QD1G.Hm[N] / (QD1G.Hp[N] + QD1G.CbfR * hxe[N]);
            }

            const double epsk1g = 1.0e-6;
            const double epsphi1g1 = 1.0e-6;
            const double epsphi1g2 = 1.0e-15;

            phi1gold = QD1G.phibar;
            int l = 0;
            while (l < core.maxin){
                l++;
                double lambda1gold = solver.lambda;

                // reset and fill coefficient matrix
                Al.assign(n, 0.0);
                Ad.assign(n, 0.0);
                Au.assign(n, 0.0);
                for (int i = 0 ; i < N ; i++){
                    Ad[i] = a2[i];
                    if (i > 0){
                        Al[i] = a1[i];
                    }
                    if (i < N - 1){
                        Au[i] = a3[i];
                    }
                }

                // fission source from previous iterate
                for (int i = 0 ; i < N ; i++){
                    invlambda = 1.0/ lambda1gold;
                    b[i] = invlambda * QD1G.nu[i] * QD1G.xsfiss[i] * mesh.hx[i];
                    bvec[i] = b[i] * phi1gold[i];
                }

                // solve linear system with the Thomas algorithm
                for (int k = 1; k < n; k++){
                    const double f = Al[k] / Ad[k-1];
                    Ad[k] -= f * Au[k-1];
                    Al[k] = 0.0;
                    bvec[k] -= f * bvec[k-1];
                }
                // back substitution
                phivec[n-1] = bvec[n-1] / Ad[n-1];
                for (int k = n - 2; k >= 0; k--){
                    phivec[k] = (bvec[k] - Au[k] * phivec[k+1]) / Ad[k];
                }

                // compute boundary leakages for eigenvalue
                const double JL = QD1G.CbfL * QD1G.Hp[0] * phivec[0] / (QD1G.Hm[0] - QD1G.CbfL * hxe[0]);
                const double JR = QD1G.CbfR * QD1G.Hm[N] * phivec[N - 1] / (QD1G.Hp[N] + QD1G.CbfR * hxe[N]);

                // compute new eigenvalue estimate
                double sfis = 0.0;
                double absr = 0.0;
                for (int i = 0 ; i < N ; i++){
                    sfis += QD1G.nu[i] * QD1G.xsfiss[i] * mesh.hx[i] * phivec[i];
                    absr += QD1G.xsabs[i] * mesh.hx[i] * phivec[i];
                }
                solver.lambda = sfis / (JR - JL + absr);

                // stopping criteria
                double dphi = 0.0;
                double phimax = 0.0;
                for (int i = 0 ; i < N ; i++){
                    dphi = std::max(dphi, std::fabs(phivec[i] - phi1gold[i]));
                    phimax = std::max(phimax, std::fabs(phivec[i]));
                }
                phi1gold = phivec;
                if (dphi <= epsphi1g1 * phimax + epsphi1g2 && std::fabs(solver.lambda - lambda1gold) <= epsk1g * solver.lambda){
                    break;
                }
            }

            // outer iteration errors on the one-group solution
            maxepsflx = 0.0;
            maxsfeps = 0.0;
            phimax = 0.0;
            for (int i = 0 ; i < N ; i++){
                maxepsflx = std::max(maxepsflx, std::fabs(phivec[i] - QD1G.phibar[i]));
                phimax = std::max(phimax, std::fabs(phivec[i]));
                sfisnew[i] = QD1G.nu[i] * QD1G.xsfiss[i] * phivec[i];
                sfiseps[i] = std::fabs(sfisnew[i] - sfisold[i]) / std::max(std::fabs(sfisold[i]), 1.0e-14);
                maxsfeps = std::max(maxsfeps, sfiseps[i]);
            }

            // store one group solution
            QD1G.phibar = phivec;
            QD1G.k = solver.lambda;

            // reconstruct edge fluxes and currents
            QD1G.phie[0] = QD1G.Hp[0] * phivec[0] / (QD1G.Hm[0] - QD1G.CbfL * hxe[0]);
            QD1G.phie[N] = QD1G.Hm[N] * phivec[N - 1] / (QD1G.Hp[N] + QD1G.CbfR * hxe[N]);
            QD1G.je[0] = QD1G.CbfL * QD1G.phie[0];
            QD1G.je[N] = QD1G.CbfR * QD1G.phie[N];
            for (int i = 1 ; i < N ; i++){
                QD1G.phie[i] = QD1G.Am[i] * phivec[i - 1] + QD1G.Ap[i] * phivec[i];
                QD1G.je[i] = (QD1G.Hm[i] * phivec[i - 1] - QD1G.Hp[i] * phivec[i]) / hxe[i];
            }

        } else {
            // compute new keff from new fission source via standard outer iteration
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

            solver.lambda = kold * (sfistotnew / sfistotold);
        }

        rhoflx = maxepsflx / maxepsflxlm2;
        keps = std::fabs(solver.lambda - kold);
        rhok = keps / kepslm2;

        solver.keff.push_back(solver.lambda);
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
        bool converged = false;
        if (core.acceleration == "MLQD"){
            const double epskout = 1.0e-5;
            const double epsphiout1 = 1.0e-5;
            const double epsphiout2 = 1.0e-15;
            converged = (maxepsflx <= epsphiout1 * phimax + epsphiout2) && (keps <= epskout * solver.lambda);
        } else {
            stopkcrit = core.epsk * (1.0 / rhok - 1.0);
            stopflx = core.epsflx * (1.0 / rhoflx - 1.0);
            converged = (keps < stopkcrit && maxepsflx < stopflx && maxsfeps < core.epspow);
        }

        if (converged) {
            if (core.search != "critical_size") {
                std::cout << "Outer iteration convergence reached in " << lout << " iterations.\n";
                if (core.acceleration != "MLQD"){
                    std::cout << "Final convergence criterion: \n";
                    std::cout << "  keff" << std::setw(16) << std::scientific << std::setprecision(6) << stopkcrit << "\n";
                    std::cout << "  flux" << std::setw(16) << std::scientific << std::setprecision(6) << stopflx << "\n";
                    std::cout << "\n";
                }
            }
            break;
        }
    }

    // normalize
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

    if (core.acceleration == "MLQD"){
        double phitotLO = 0.0;
        for (int g = 0 ; g < ngroups ; g++){
            for (int i = 0 ; i < mesh.ncells ; i++){
                phitotLO += QDMG.phibar[g][i] * mesh.hx[i];
            }
        }
        const double normLO = H / phitotLO;
        for (int g = 0 ; g < ngroups ; g++){
            for (int i = 0 ; i < mesh.ncells ; i++){
                QDMG.phibar[g][i] *= normLO;
            }
            for (int i = 0 ; i < mesh.nedges ; i++){
                QDMG.phie[g][i] *= normLO;
                QDMG.je[g][i] *= normLO;
            }
        }
        for (int i = 0 ; i < mesh.ncells ; i++){
            QD1G.phibar[i] *= normLO;
        }
        for (int i = 0 ; i < mesh.nedges ; i++){
            QD1G.phie[i] *= normLO;
            QD1G.je[i] *= normLO;
        }
    }

    solver.k_outer += lout;
}