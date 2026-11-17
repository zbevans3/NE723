#include "source_iteration.h"
#include "transport_sweep.h"

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <string>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <iomanip>


void SourceIteration::solve(
        coreData& core,
        quadratureData& quadrature,
        std::vector<materialData>& materials,
        meshData& mesh,
        solverData& solver,
        QDData& QD,
        QDMGData& QDMG,
        QD1GData& QD1G,
        double invlambda
    )    
{
    int l = 0;
    double maxeps = 1.0e6;
    double maxepslm2 = 1.0e6;
    double rhol = 1.0e6;
    double maxsres = 1.0e6;
    double maxsres1 = 1.0e6;
    double maxires = 1.0e6;
    double maxires1 = 1.0e6;
    double stopcrit = 1.0e6;

    phibarold = solver.phibar;
    phihatold = solver.phihat;
    const int ngroups = core.ngroups;
    Qi.assign(ngroups, std::vector<double>(mesh.ncells, 0.0));
    Qhati.assign(ngroups, std::vector<double>(mesh.ncells, 0.0));
    phieps.assign(ngroups, std::vector<double>(mesh.ncells, 1.0e6));
    qfis.assign(ngroups, std::vector<double>(mesh.ncells, 0.0));
    qfishat.assign(ngroups, std::vector<double>(mesh.ncells, 0.0));

    TransportSweep sweeper;

    if (core.acceleration == "QD"){
      rhol = 0.5;
    }

    // assign coefficient that determines LD or LLD
    double alpha = (core.scheme == "LD") ? 1.0 : 3.0;
    if (core.search == "fixed_source" && core.itdebug){
      std::cout << std::right << std::setw(12) << "  SI" << std::setw(16) << "Flux error" << std::setw(16) << "S flux res"
                << std::setw(16) << "S moment res" << std::setw(16) << "I flux res"
                << std::setw(16) << "I moment res" << "\n";
    };

    if (core.acceleration != "MLQD"){
        // compute initial fission source
        // if eigenvalue calculation, then fission source is updated only on outer iteration
        for (int i = 0 ; i < mesh.ncells ; i++){
            // total fission rate in cell i, summed over all groups
            double fis = 0.0;
            double fishat = 0.0;
            for (int gp = 0 ; gp < ngroups ; gp++){
                fis += materials[mesh.mapmat[i]].nu[gp] * materials[mesh.mapmat[i]].xsfiss[gp] * solver.phibar[gp][i];
                fishat += materials[mesh.mapmat[i]].nu[gp] * materials[mesh.mapmat[i]].xsfiss[gp] * solver.phihat[gp][i];
            }
            // compute group contribution based on chi
            for (int g = 0 ; g < ngroups ; g++){
                qfis[g][i] = invlambda * materials[mesh.mapmat[i]].xschi[g] * fis;
                qfishat[g][i] = invlambda * materials[mesh.mapmat[i]].xschi[g] * fishat;
            }
        }

        while (l < core.maxin){
            l++;
            // update l-2 values
            (l == 1) ? maxepslm2 = 1.0 : maxepslm2 = maxeps;

            // update l-1 values
            phibarold = solver.phibar;
            phihatold = solver.phihat;

            // fission source is updated every inner iteration only for fixed source problems
            if (core.search == "fixed_source"){
                for (int i = 0 ; i < mesh.ncells ; i++){
                    double fis = 0.0;
                    double fishat = 0.0;
                    for (int gp = 0 ; gp < ngroups ; gp++){
                        fis += materials[mesh.mapmat[i]].nu[gp] * materials[mesh.mapmat[i]].xsfiss[gp] * solver.phibar[gp][i];
                        fishat += materials[mesh.mapmat[i]].nu[gp] * materials[mesh.mapmat[i]].xsfiss[gp] * solver.phihat[gp][i];
                    }
                    for (int g = 0 ; g < ngroups ; g++){
                        qfis[g][i] = invlambda * materials[mesh.mapmat[i]].xschi[g] * fis;
                        qfishat[g][i] = invlambda * materials[mesh.mapmat[i]].xschi[g] * fishat;
                    }
                }
            }

            // sweep group g, then update its fluxes
            for (int g = 0 ; g < ngroups ; g++){
              // computed l-th iterate source and 1st moment for group g
              for (int i = 0 ; i < mesh.ncells ; i++){
                  double qexternal = core.search == "fixed_source" ? materials[mesh.mapmat[i]].q[g] : 0.0;
                  // scattering into group g
                  double scat = 0.0;
                  double scathat = 0.0;
                  for (int gp = 0 ; gp < ngroups ; gp++){
                      scat += materials[mesh.mapmat[i]].xsscatmat[gp][g] * solver.phibar[gp][i];
                      scathat += materials[mesh.mapmat[i]].xsscatmat[gp][g] * solver.phihat[gp][i];
                  }
                  Qi[g][i] = 0.5 * mesh.hx[i] * (scat + qfis[g][i] + qexternal);
                  Qhati[g][i] = 0.5 * mesh.hx[i] * (scathat + qfishat[g][i]);
              }

              // perform mesh sweep for each ordinate
              if ((l > 1 && core.acceleration == "QD") || core.acceleration == "none"){
                sweeper.sweep_group(g, Qi[g], Qhati[g], alpha, quadrature, materials, mesh, solver);
              }

              // update group g fluxes
              update_flux_and_current(
                      quadrature,
                      mesh,
                      solver,
                      g
              );
            }

            // compute residuals of the transport solution
            if (core.acceleration == "none" || l > 1){
              compute_residuals(core, materials, mesh, solver, alpha, invlambda);
            }

            // if (core.acceleration == "QD"){
            //   // do QD acceleration
            //   QD_acceleration(
            //           core,
            //           quadrature,
            //           materials,
            //           mesh,
            //           solver,
            //           QD,
            //           invlambda,
            //           l > 1
            //   );
            //   // do prolongation on scalar flux and first legendre moment
            //   for (int i = 0 ; i < mesh.ncells ; i++){
            //     solver.phibar[i] = QD.phibarLO[i];
            //     solver.phihat[i] = (QD.phieLO[i+1] - QD.phieLO[i]) * 0.5;
            //     phieps[i] = fabs(solver.phibar[i]-phibarold[i]);
            //   }
            // };

            // max norms over all groups and cells
            maxeps = 0.0;
            maxsres = 0.0;
            maxsres1 = 0.0;
            maxires = 0.0;
            maxires1 = 0.0;
            for (int g = 0 ; g < ngroups ; g++){
                for (int i = 0 ; i < mesh.ncells ; i++){
                    maxeps = std::max(maxeps, fabs(phieps[g][i]));
                    maxsres = std::max(maxsres, fabs(solver.sresidual[g][i]));
                    maxsres1 = std::max(maxsres1, fabs(solver.sresidual1[g][i]));
                    maxires = std::max(maxires, fabs(solver.iresidual[g][i]));
                    maxires1 = std::max(maxires1, fabs(solver.iresidual1[g][i]));
                }
            }
            if (l > 2){
                rhol = maxeps / maxepslm2;
            } else {
                rhol = 0.5;
            }

            if (core.itdebug){
              std::cout << std::right << std::scientific << std::setprecision(6)
                        << std::setw(12) << l
                        << std::setw(16) << maxeps
                        << std::setw(16) << maxsres
                        << std::setw(16) << maxsres1
                        << std::setw(16) << maxires
                        << std::setw(16) << maxires1
                        << '\n';
            }

            solver.ilinfphi.push_back(maxeps);
            solver.ispecrad.push_back(rhol);

            stopcrit = core.epsflx * (1.0 / rhol - 1.0);
            if (maxeps < stopcrit) {
                if (core.search != "eigenvalue" && core.search != "critical_size"){
                  std::cout << "Source iteration convergence reached in " << l << " iterations.\n";
                  std::cout << "Final convergence criterion: "
                            << std::scientific << std::setprecision(6) << stopcrit << "\n";
                }
                break;
            }
        }

        if (core.acceleration == "QD"){
            // compute measure of transport effects
            for (int i = 0; i < mesh.ncells ; i++){
                const double invhx = 1.0 / mesh.hx[i];
                QD.mLO[i] = fabs(1.0/3.0 - QD.Eddbar[i]);
                QD.mprimeLO[i] = fabs((QD.Edde[i+1]-QD.Edde[i]) * invhx);
            }
        }
        solver.k_inner += l;
    } else {
        // for MLQD, there are no source iterations
        // build fission source from one group scalar flux
        for (int i = 0 ; i < mesh.ncells ; i++){
            // total fission rate in cell i, summed over all groups
            double fis = 0.0;
            double fishat = 0.0;
            for (int gp = 0 ; gp < ngroups ; gp++){
                fis += materials[mesh.mapmat[i]].nu[gp] * materials[mesh.mapmat[i]].xsfiss[gp] * solver.phibar[gp][i];
                fishat += materials[mesh.mapmat[i]].nu[gp] * materials[mesh.mapmat[i]].xsfiss[gp] * solver.phihat[gp][i];
            }
            // compute group contribution based on chi
            for (int g = 0 ; g < ngroups ; g++){
                qfis[g][i] = invlambda * materials[mesh.mapmat[i]].xschi[g] * fis;
                qfishat[g][i] = invlambda * materials[mesh.mapmat[i]].xschi[g] * fishat;
            }
        }

    }
}

void SourceIteration::update_flux_and_current(
        quadratureData& quadrature,
        meshData& mesh,
        solverData& solver,
        int g
){

    // compute new cell edge scalar flux and current for group g
    for (int i = 0 ; i < mesh.nedges ; i++){
       solver.phie[g][i] = 0.0;
       solver.je[g][i] = 0.0;
       for (int u = 0 ; u < quadrature.order ; u++){
         solver.phie[g][i] += quadrature.weight[u] * solver.psie[g][u][i];
         solver.je[g][i] += quadrature.angle[u] * quadrature.weight[u] * solver.psie[g][u][i];
       }
    }

    // compute new cell avg scalar flux and first moment scalar flux for group g
    for (int i = 0 ; i < mesh.ncells ; i++){
       solver.phibar[g][i] = 0.0;
       solver.phihat[g][i] = 0.0;
       solver.jbar[g][i] = 0.0;
       solver.jhat[g][i] = 0.0;

       // compute new flux and currents
       for (int u = 0 ; u < quadrature.order ; u++){
         solver.phibar[g][i] += quadrature.weight[u] * solver.psibar[g][u][i];
         solver.phihat[g][i] += quadrature.weight[u] * solver.psihat[g][u][i];
         solver.jbar[g][i] += quadrature.angle[u] * quadrature.weight[u] * solver.psibar[g][u][i];
         solver.jhat[g][i] += quadrature.angle[u] * quadrature.weight[u] * solver.psihat[g][u][i];
       }
       phieps[g][i] = fabs(solver.phibar[g][i]-phibarold[g][i]);
    }

}

void SourceIteration::compute_residuals(
        coreData& core,
        std::vector<materialData>& materials,
        meshData& mesh,
        solverData& solver,
        double& alpha,
        double invlambda
){
    const int ngroups = core.ngroups;

    for (int i = 0 ; i < mesh.ncells ; i++){

       // total fission rate in cell i from the current fluxes (fixed source only;
       // for eigenvalue problems the fission source is held fixed during inner iterations)
       double fis = 0.0;
       double fishat = 0.0;
       if (core.search == "fixed_source"){
           for (int gp = 0 ; gp < ngroups ; gp++){
               fis += materials[mesh.mapmat[i]].nu[gp] * materials[mesh.mapmat[i]].xsfiss[gp] * solver.phibar[gp][i];
               fishat += materials[mesh.mapmat[i]].nu[gp] * materials[mesh.mapmat[i]].xsfiss[gp] * solver.phihat[gp][i];
           }
       }

       for (int g = 0 ; g < ngroups ; g++){
          const double sigt = materials[mesh.mapmat[i]].xstot[g];
          const double hx = mesh.hx[i];

          // scheme: balance against the source used in this sweep (Qi = 0.5*hx*source)
          solver.sresidual[g][i] = solver.je[g][i+1] - solver.je[g][i] + sigt * hx * solver.phibar[g][i]
                                   - 2.0 * Qi[g][i];
          solver.sresidual1[g][i] = 3.0*(solver.je[g][i+1] + solver.je[g][i] - 2.0*solver.jbar[g][i])
                                    + alpha * sigt * hx * solver.phihat[g][i]
                                    - alpha * 2.0 * Qhati[g][i];

          // iterative: balance against the source rebuilt from the current fluxes
          double scat = 0.0;
          double scathat = 0.0;
          for (int gp = 0 ; gp < ngroups ; gp++){
              scat += materials[mesh.mapmat[i]].xsscatmat[gp][g] * solver.phibar[gp][i];
              scathat += materials[mesh.mapmat[i]].xsscatmat[gp][g] * solver.phihat[gp][i];
          }
          double qf = qfis[g][i];
          double qfhat = qfishat[g][i];
          double qexternal = 0.0;
          if (core.search == "fixed_source"){
              qf = invlambda * materials[mesh.mapmat[i]].xschi[g] * fis;
              qfhat = invlambda * materials[mesh.mapmat[i]].xschi[g] * fishat;
              qexternal = materials[mesh.mapmat[i]].q[g];
          }

          solver.iresidual[g][i] = solver.je[g][i+1] - solver.je[g][i] + sigt * hx * solver.phibar[g][i]
                                   - (scat + qf + qexternal) * hx;
          solver.iresidual1[g][i] = 3.0*(solver.je[g][i+1] + solver.je[g][i] - 2.0*solver.jbar[g][i])
                                    + alpha * sigt * hx * solver.phihat[g][i]
                                    - alpha * (scathat + qfhat) * hx;
       }
    }

}

// void SourceIteration::QD_acceleration(
//         coreData& core,
//         quadratureData& quadrature,
//         std::vector<materialData>& materials,
//         meshData& mesh,
//         solverData& solver,
//         QDData& QD,
//         double invlambda,
//         bool update_factors
// ){

//   // compute E, C_L, C_R from the high-order solution
//   if (update_factors){
//     // boundary factors
//     double nsumL = 0.0;
//     double dsumL = 0.0;
//     double nsumR = 0.0;
//     double dsumR = 0.0;
//     for (int u = 0; u < quadrature.order; u++){
//         const double weight  = quadrature.weight[u];
//         const double mu = quadrature.angle[u];
//         if (mu < 0.0){
//             // leaving through the left boundary
//             nsumL += weight * mu * solver.psie[u].front();
//             dsumL += weight * solver.psie[u].front();
//         } else {
//             // leaving through the right boundary
//             nsumR += weight * mu * solver.psie[u].back();
//             dsumR += weight * solver.psie[u].back();
//         }
//     }
//     QD.CbfL = nsumL / dsumL;
//     QD.CbfR = nsumR / dsumR;

//     // Eddington factors at cell edges
//     for (int i = 0; i < mesh.nedges; i++){
//         double nsum = 0.0;
//         double dsum = 0.0;
//         for (int u = 0; u < quadrature.order; u++){
//             const double mu = quadrature.angle[u];
//             const double weight = quadrature.weight[u];
//             nsum += weight * mu * mu * solver.psie[u][i];
//             dsum += weight * solver.psie[u][i];
//         }
//         QD.Edde[i] = nsum / dsum;
//     }

//     // Eddington factors from cell averages
//     for (int i = 0; i < mesh.ncells; i++){
//         double nsum = 0.0;
//         double dsum = 0.0;
//         for (int u = 0; u < quadrature.order; u++){
//             const double mu = quadrature.angle[u];
//             const double weight = quadrature.weight[u];
//             nsum += weight * mu * mu * solver.psibar[u][i];
//             dsum += weight * solver.psibar[u][i];
//         }
//         QD.Eddbar[i] = nsum / dsum;
//     }
//   }

//     // build linear system
//     const int N = mesh.ncells;
//     const int n = 2 * N + 1;

//     std::vector<std::vector<double>> Amat(n, std::vector<double>(n, 0.0));
//     std::vector<double> bvec(n, 0.0);
//     std::vector<double> phivec(n, 0.0);

//     // computencoefficients
//     std::vector<double> a1(N);
//     std::vector<double> a2(N);
//     std::vector<double> b(N);
//     for (int i = 0; i < N; i++){
//         double sigmar = materials[mesh.mapmat[i]].xstot - materials[mesh.mapmat[i]].xsscat;
//         double q = 0.0;
//         if (core.search == "fixed_source"){
//             sigmar -= invlambda * materials[mesh.mapmat[i]].nu * materials[mesh.mapmat[i]].xsfiss;
//             q = materials[mesh.mapmat[i]].q;
//         } else {
//             q = qfis[i];
//         }
//         a1[i] = 2.0 / (materials[mesh.mapmat[i]].xstot * mesh.hx[i]);
//         a2[i] = 2.0 * a1[i] * QD.Eddbar[i] + sigmar * mesh.hx[i];
//         b[i] = q * mesh.hx[i];
//     }

//     // compute incoming flux and current, summed over groups (total for the one-group LO problem)
//     double phiL = 0.0;
//     double JinL = 0.0;
//     double phiR = 0.0;
//     double JinR = 0.0;
//     for (int g = 0; g < core.ngroups; g++){
//         for (int u = 0; u < quadrature.order; u++){
//             const double weight = quadrature.weight[u];
//             const double mu = quadrature.angle[u];
//             if (mu > 0.0){
//                 phiL += weight * psilbc[g][u];
//                 JinL += weight * mu * psilbc[g][u];
//             } else {
//                 phiR += weight * psirbc[g][u];
//                 JinR += weight * mu * psirbc[g][u];
//             }
//         }
//     }

//     // left edge
//     Amat[0][0] = a1[0] * QD.Edde[0];
//     Amat[0][1] = -a1[0] * QD.Eddbar[0];
//     if (mesh.mapbc.front() != 2){
//         Amat[0][0] -= QD.CbfL;
//         bvec[0] = JinL - QD.CbfL * phiL;
//     }

//     // balance equation
//     for (int i = 0; i < N; i++){
//         const int k = 2 * i + 1;
//         Amat[k][k-1] = -a1[i] * QD.Edde[i];
//         Amat[k][k] = a2[i];
//         Amat[k][k+1] = -a1[i] * QD.Edde[i+1];
//         bvec[k] = b[i];
//     }

//     // current continuity on interior
//     for (int i = 0; i < N - 1; i++){
//         const int k = 2 * i + 2;
//         Amat[k][k-1] = a1[i] * QD.Eddbar[i];
//         Amat[k][k] = -(a1[i] + a1[i+1]) * QD.Edde[i+1];
//         Amat[k][k+1] = a1[i+1] * QD.Eddbar[i+1];
//     }

//     // right edge
//     Amat[n-1][n-2] = a1[N-1] * QD.Eddbar[N-1];
//     Amat[n-1][n-1] = -a1[N-1] * QD.Edde[N];
//     if (mesh.mapbc.back() != 2){
//         Amat[n-1][n-1] -= QD.CbfR;
//         bvec[n-1] = JinR - QD.CbfR * phiR;
//     }

//     // solve linear system with the Thomas algorithm
//     for (int k = 1; k < n; k++){
//         const double f = Amat[k][k-1] / Amat[k-1][k-1];
//         Amat[k][k] -= f * Amat[k-1][k];
//         Amat[k][k-1] = 0.0;
//         bvec[k] -= f * bvec[k-1];
//     }
//     // back substitution
//     phivec[n-1] = bvec[n-1] / Amat[n-1][n-1];
//     for (int k = n - 2; k >= 0; k--){
//         phivec[k] = (bvec[k] - Amat[k][k+1] * phivec[k+1]) / Amat[k][k];
//     }

//     // update LO scalar flux values
//     for (int i = 0; i < N; i++) {
//         QD.phieLO[i]   = phivec[2 * i];
//         QD.phibarLO[i] = phivec[2 * i + 1];
//     }
//     QD.phieLO[N] = phivec[2 * N];

//     // update LO current values
//     for (int i = 0; i < N; i++) {
//         QD.jeLO[i] = a1[i] * (QD.Edde[i] * QD.phieLO[i] - QD.Eddbar[i] * QD.phibarLO[i]);
//     }
//     QD.jeLO[N] = a1[N-1] * (QD.Eddbar[N-1] * QD.phibarLO[N-1] - QD.Edde[N] * QD.phieLO[N]);

// }