#include "source_iteration.h"

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
        double invlambda
    )    
{

    phibarold = solver.phibar;
    phihatold = solver.phihat;
    Qi.assign(mesh.ncells, 0.0);
    Qhati.assign(mesh.ncells, 0.0);
    psilbc.assign(quadrature.order, 0.0);
    psirbc.assign(quadrature.order, 0.0);
    phieps.assign(mesh.ncells, 1.0e6);

    // fill left boundary flux array
    if (mesh.mapbc.front() == 1){
      // vacuum bc
      psilbc.assign(quadrature.order, 0.0);
    } else if (mesh.mapbc.front() == 2) {
      // reflective bc, initalize to 0
      psilbc.assign(quadrature.order, 0.0);
    } else if (mesh.mapbc.front() == 3){
      // incoming
      psilbc.assign(quadrature.order, core.psileft);
    }
    
    // fill right boundary flux array
    if (mesh.mapbc.back() == 1){
      // vacuum bc
      psirbc.assign(quadrature.order, 0.0);
    } else if (mesh.mapbc.back() == 2) {
      // reflective bc, initialize to 0
      psirbc.assign(quadrature.order, 0.0);
    } else if (mesh.mapbc.back() == 3){
      // incoming
      psirbc.assign(quadrature.order, core.psiright);
    }
    if (core.search == "fixed_source" && core.itdebug){
      std::cout << std::right << std::setw(12) << "  SI" << std::setw(16) << "Flux error" << std::setw(16) << "S flux res"
                << std::setw(16) << "S moment res" << std::setw(16) << "I flux res"
                << std::setw(16) << "I moment res" << "\n";
    };


    // compute initial fission source
    // if eigenvalue calculation, then fission source is updated only on outer iteration
    std::vector<double> qfis(mesh.ncells, 0.0);
    std::vector<double> qfishat(mesh.ncells, 0.0);
    
    for (int i = 0; i < mesh.ncells; ++i) {  
        qfis[i] = invlambda * materials[mesh.mapmat[i]].nu *
                  materials[mesh.mapmat[i]].xsfiss * solver.phibar[i];
        qfishat[i] = invlambda * materials[mesh.mapmat[i]].nu *
                     materials[mesh.mapmat[i]].xsfiss * solver.phihat[i];
    }

    int l = 0;
    double betai = 0;
    double Gi = 0.0;
    double invmu = 0.0;
    double invGi = 0.0;
    double maxeps = 1.0e6;
    double maxepslm2 = 1.0e6;
    double rhol = 1.0e6;
    double maxsres = 1.0e6;
    double maxsres1 = 1.0e6;
    double maxires = 1.0e6;
    double maxires1 = 1.0e6;
    double stopcrit = 1.0e6;
    while (l < core.maxin){
        l++;
        // update l-2 values
        (l == 1) ? maxepslm2 = 1.0 : maxepslm2 = maxeps;

        // update l-1 values
        phibarold = solver.phibar;
        phihatold = solver.phihat;

        // computed l-th iterate source and 1st moment for cell i
        for (int i = 0 ; i < mesh.ncells ; i++){
          if (core.search == "fixed_source"){
            qfis[i] = invlambda * materials[mesh.mapmat[i]].nu *
                      materials[mesh.mapmat[i]].xsfiss * solver.phibar[i];
            qfishat[i] = invlambda * materials[mesh.mapmat[i]].nu *
                         materials[mesh.mapmat[i]].xsfiss * solver.phihat[i];
          }
          double qexternal = core.search == "fixed_source" ? materials[mesh.mapmat[i]].q : 0.0;
          
          Qi[i] = 0.5 * mesh.hx[i] * (materials[mesh.mapmat[i]].xsscat * solver.phibar[i] 
                                      + qfis[i] + qexternal);
          Qhati[i] = 0.5 * mesh.hx[i] * (materials[mesh.mapmat[i]].xsscat * solver.phihat[i] 
                                      + qfishat[i]);
        }

        // perform mesh sweep for each ordinate
        for (int u = 0 ; u < quadrature.order ; u++){
           invmu = 1.0 / quadrature.angle[u];
           if (quadrature.angle[u] > 0){
              // positive sweep (x+)
              solver.psie[u].front() = psilbc[u];
              
              // step through cells
              for (int i = 0 ; i < mesh.ncells ; i++){                                                    
                betai = materials[mesh.mapmat[i]].xstot * mesh.hx[i] * invmu;
                Gi = betai*betai + 4.0 * betai + 6.0;
                invGi = 1.0/Gi;

                // compute cell  avg angular flux and first legendre moment
                solver.psibar[u][i] =  ((3.0 + betai) * Qi[i] - Qhati[i]) * (invmu * invGi)
                                         + (6.0 + betai) * invGi * solver.psie[u][i];
                solver.psihat[u][i] =  (3.0 * Qi[i] + (1.0 + betai) * Qhati[i]) * (invmu * invGi)
                                         - (3.0 * betai) * invGi * solver.psie[u][i];

                // compute outgoing angular flux
                solver.psie[u][i+1] = solver.psibar[u][i] + solver.psihat[u][i];                 
              }
           }
           else{
              // negative sweep (x-)
              solver.psie[u].back() = psirbc[u];

              // step through cells
              for (int i = mesh.ncells - 1; i >= 0; --i) {                                                 
                betai = materials[mesh.mapmat[i]].xstot * mesh.hx[i] * invmu;
                Gi = betai*betai - 4.0 * betai + 6.0;
                invGi = 1.0/Gi;

                // compute cell  avg angular flux and first legendre moment
                solver.psibar[u][i] =  ((betai - 3.0) * Qi[i] - Qhati[i]) * (invmu * invGi)
                                         + (6.0 - betai) * invGi * solver.psie[u][i+1];
                solver.psihat[u][i] =  (3.0 * Qi[i] + (betai - 1.0) * Qhati[i]) * (invmu * invGi)
                                         - (3.0 * betai) * invGi * solver.psie[u][i+1];

                // compute outgoing angular flux
                solver.psie[u][i] = solver.psibar[u][i] - solver.psihat[u][i];                 
              }
           }
        }

        // update reflective edge fluxes with opposite ordinate's
        if (mesh.mapbc.front() == 2 || mesh.mapbc.back() == 2){
          for (int u = 0; u < quadrature.order; ++u) {
              if (mesh.mapbc.front() == 2 && quadrature.angle[u] > 0.0) {
                   psilbc[u] = solver.psie[quadrature.reflmap[u]].front();
              }
              if (mesh.mapbc.back() == 2 && quadrature.angle[u] < 0.0) {
                  psirbc[u] = solver.psie[quadrature.reflmap[u]].back();
              }
          }
        }

        // compute new cell edge scalar flux and current
        for (int i = 0 ; i < mesh.nedges ; i++){
           solver.phie[i] = 0.0;
           solver.je[i] = 0.0;
           for (int u = 0 ; u < quadrature.order ; u++){
             solver.phie[i] += quadrature.weight[u] * solver.psie[u][i];
             solver.je[i] += quadrature.angle[u] * quadrature.weight[u] * solver.psie[u][i];
           }
        }

        // compute new cell avg scalar flux and first moment scalar flux
        for (int i = 0 ; i < mesh.ncells ; i++){
           solver.phibar[i] = 0.0;
           solver.phihat[i] = 0.0;
           solver.jbar[i] = 0.0;
           solver.jhat[i] = 0.0;
                             
           // compute new flux and currents
           for (int u = 0 ; u < quadrature.order ; u++){
             solver.phibar[i] += quadrature.weight[u] * solver.psibar[u][i];
             solver.phihat[i] += quadrature.weight[u] * solver.psihat[u][i];
             solver.jbar[i] += quadrature.angle[u] * quadrature.weight[u] * solver.psibar[u][i];
             solver.jhat[i] += quadrature.angle[u] * quadrature.weight[u] * solver.psihat[u][i];
           }              
           phieps[i] = fabs(solver.phibar[i]-phibarold[i]);
        }
        
        // compute residuals
        for (int i = 0 ; i < mesh.ncells ; i++){
           // scheme
           solver.sresidual[i] = solver.je[i+1] - solver.je[i] + materials[mesh.mapmat[i]].xstot * mesh.hx[i] * 
                                 solver.phibar[i] - (materials[mesh.mapmat[i]].xsscat + materials[mesh.mapmat[i]].nu *
                                 materials[mesh.mapmat[i]].xsfiss) * mesh.hx[i] * phibarold[i] - 
                                 materials[mesh.mapmat[i]].q * mesh.hx[i];
           solver.sresidual1[i] = 3.0*(solver.je[i+1] + solver.je[i] - 2.0*solver.jbar[i]) +
                                  materials[mesh.mapmat[i]].xstot * mesh.hx[i] * solver.phihat[i] - 
                                  (materials[mesh.mapmat[i]].xsscat + materials[mesh.mapmat[i]].nu * materials[mesh.mapmat[i]].xsfiss) *
                                  mesh.hx[i] * phihatold[i];
           // iterative
           solver.iresidual[i] = solver.je[i+1] - solver.je[i] + materials[mesh.mapmat[i]].xstot * mesh.hx[i] * 
                                 solver.phibar[i] - (materials[mesh.mapmat[i]].xsscat + materials[mesh.mapmat[i]].nu *
                                 materials[mesh.mapmat[i]].xsfiss) * mesh.hx[i] * solver.phibar[i] - 
                                 materials[mesh.mapmat[i]].q * mesh.hx[i];
           solver.iresidual1[i] = 3.0*(solver.je[i+1] + solver.je[i] - 2.0*solver.jbar[i]) +
                                  materials[mesh.mapmat[i]].xstot * mesh.hx[i] * solver.phihat[i] - 
                                  (materials[mesh.mapmat[i]].xsscat + materials[mesh.mapmat[i]].nu * materials[mesh.mapmat[i]].xsfiss) *
                                  mesh.hx[i] * solver.phihat[i];
           solver.sresidual[i] = std::abs(solver.sresidual[i]);
           solver.sresidual1[i] = std::abs(solver.sresidual1[i]);
           solver.iresidual[i] = std::abs(solver.iresidual[i]);
           solver.iresidual1[i] = std::abs(solver.iresidual1[i]);                                  
        }

        maxeps = *std::max_element(phieps.begin(), phieps.end());
        rhol = maxeps / maxepslm2;
        maxsres = *std::max_element(solver.sresidual.begin(), solver.sresidual.end());
        maxsres1 = *std::max_element(solver.sresidual1.begin(), solver.sresidual1.end());
        maxires = *std::max_element(solver.iresidual.begin(), solver.iresidual.end());
        maxires1 = *std::max_element(solver.iresidual1.begin(), solver.iresidual1.end());

        if (core.itdebug){
          std::cout << std::right << std::scientific << std::setprecision(6)
                    << std::setw(12)  << l
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

    solver.k_inner += l;
}