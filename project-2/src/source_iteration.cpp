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
        QDData& QD,
        double invlambda
    )    
{
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

    phibarold = solver.phibar;
    phihatold = solver.phihat;
    Qi.assign(mesh.ncells, 0.0);
    Qhati.assign(mesh.ncells, 0.0);
    psilbc.assign(quadrature.order, 0.0);
    psirbc.assign(quadrature.order, 0.0);
    phieps.assign(mesh.ncells, 1.0e6);

    if (core.acceleration == "QD"){
      // build QD data
      QD.CbfL = -0.5;
      QD.CbfR = 0.5;
      QD.Edde.assign(mesh.nedges, 1.0/3.0);
      QD.Eddbar.assign(mesh.ncells, 1.0/3.0);
      rhol = 0.5;
      QD.phibarLO.assign(mesh.ncells, 0.0);
      QD.phieLO.assign(mesh.nedges, 0.0);
      QD.jeLO.assign(mesh.nedges, 0.0);
      QD.mLO.assign(mesh.ncells, 0.0);
      QD.mprimeLO.assign(mesh.ncells, 0.0);
    }

    // assign coefficient that determines LD or LLD
    double alpha = (core.scheme == "LD") ? 1.0 : 3.0;

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
    qfis.assign(mesh.ncells, 0.0);
    qfishat.assign(mesh.ncells, 0.0);

    for (int i = 0; i < mesh.ncells; ++i) {  
        qfis[i] = invlambda * materials[mesh.mapmat[i]].nu *
                  materials[mesh.mapmat[i]].xsfiss * solver.phibar[i];
        qfishat[i] = invlambda * materials[mesh.mapmat[i]].nu *
                     materials[mesh.mapmat[i]].xsfiss * solver.phihat[i];
    }

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
        if ((l > 1 && core.acceleration == "QD") || core.acceleration == "none"){
            for (int u = 0 ; u < quadrature.order ; u++){
               invmu = 1.0 / quadrature.angle[u];
               if (quadrature.angle[u] > 0){
                  // positive sweep (x+)
                  solver.psie[u].front() = psilbc[u];
            
                  // step through cells
                  for (int i = 0 ; i < mesh.ncells ; i++){                                                    
                    betai = materials[mesh.mapmat[i]].xstot * mesh.hx[i] * invmu;
                    Gi = alpha*betai*betai + (3.0+alpha)*betai + 6.0;
                    invGi = 1.0/Gi;
                
                    // compute cell  avg angular flux and first legendre moment
                    solver.psibar[u][i] =  ((3.0 + alpha*betai) * Qi[i] - alpha*Qhati[i]) * (invmu * invGi)
                                             + (6.0 + alpha*betai) * invGi * solver.psie[u][i];
                    solver.psihat[u][i] =  (3.0 * Qi[i] + alpha*(1.0 + betai) * Qhati[i]) * (invmu * invGi)
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
                    Gi = alpha*betai*betai - (3.0 + alpha) * betai + 6.0;
                    invGi = 1.0/Gi;
                
                    // compute cell  avg angular flux and first legendre moment
                    solver.psibar[u][i] =  ((alpha*betai - 3.0) * Qi[i] - alpha*Qhati[i]) * (invmu * invGi)
                                             + (6.0 - alpha*betai) * invGi * solver.psie[u][i+1];
                    solver.psihat[u][i] =  (3.0 * Qi[i] + alpha*(betai - 1.0) * Qhati[i]) * (invmu * invGi)
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
        }

        update_flux_and_current(
                quadrature,
                mesh,
                solver
        );

        // compute residuals of the transport solution
        if (core.acceleration == "none" || l > 1){
          compute_residuals(materials, mesh, solver, alpha, phihatold, phibarold);
        }

        if (core.acceleration == "QD"){
          // do QD acceleration
          QD_acceleration(
                  core,
                  quadrature,
                  materials,
                  mesh,
                  solver,
                  QD,
                  invlambda,
                  l > 1
          );
          // do prolongation on scalar flux and first legendre moment
          for (int i = 0 ; i < mesh.ncells ; i++){
            solver.phibar[i] = QD.phibarLO[i];
            solver.phihat[i] = (QD.phieLO[i+1] - QD.phieLO[i]) * 0.5;
            phieps[i] = fabs(solver.phibar[i]-phibarold[i]);
          }
        };

        maxeps = *std::max_element(phieps.begin(), phieps.end());
        if (l > 2){
            rhol = maxeps / maxepslm2;
        } else {
            rhol = 0.5;
        }
        maxsres = *std::max_element(solver.sresidual.begin(), solver.sresidual.end());
        maxsres1 = *std::max_element(solver.sresidual1.begin(), solver.sresidual1.end());
        maxires = *std::max_element(solver.iresidual.begin(), solver.iresidual.end());
        maxires1 = *std::max_element(solver.iresidual1.begin(), solver.iresidual1.end());

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
}

void SourceIteration::update_flux_and_current(
        quadratureData& quadrature,
        meshData& mesh,
        solverData& solver
){

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

}

void SourceIteration::compute_residuals(
        std::vector<materialData>& materials,
        meshData& mesh,
        solverData& solver,
        double& alpha,
        std::vector<double>& phihatold,
        std::vector<double>& phibarold
){

    for (int i = 0 ; i < mesh.ncells ; i++){
       // scheme
       solver.sresidual[i] = solver.je[i+1] - solver.je[i] + materials[mesh.mapmat[i]].xstot * mesh.hx[i] * 
                             solver.phibar[i] - (materials[mesh.mapmat[i]].xsscat + materials[mesh.mapmat[i]].nu *
                             materials[mesh.mapmat[i]].xsfiss) * mesh.hx[i] * phibarold[i] - 
                             materials[mesh.mapmat[i]].q * mesh.hx[i];
       solver.sresidual1[i] = 3.0*(solver.je[i+1] + solver.je[i] - 2.0*solver.jbar[i]) +
                              alpha * materials[mesh.mapmat[i]].xstot * mesh.hx[i] * solver.phihat[i] - 
                              alpha * (materials[mesh.mapmat[i]].xsscat + materials[mesh.mapmat[i]].nu * materials[mesh.mapmat[i]].xsfiss) *
                              mesh.hx[i] * phihatold[i];
       // iterative
       solver.iresidual[i] = solver.je[i+1] - solver.je[i] + materials[mesh.mapmat[i]].xstot * mesh.hx[i] * 
                             solver.phibar[i] - (materials[mesh.mapmat[i]].xsscat + materials[mesh.mapmat[i]].nu *
                             materials[mesh.mapmat[i]].xsfiss) * mesh.hx[i] * solver.phibar[i] - 
                             materials[mesh.mapmat[i]].q * mesh.hx[i];
       solver.iresidual1[i] = 3.0*(solver.je[i+1] + solver.je[i] - 2.0*solver.jbar[i]) +
                              alpha * materials[mesh.mapmat[i]].xstot * mesh.hx[i] * solver.phihat[i] - 
                              alpha * (materials[mesh.mapmat[i]].xsscat + materials[mesh.mapmat[i]].nu * materials[mesh.mapmat[i]].xsfiss) *
                              mesh.hx[i] * solver.phihat[i];                              
    }

}

void SourceIteration::QD_acceleration(
        coreData& core,
        quadratureData& quadrature,
        std::vector<materialData>& materials,
        meshData& mesh,
        solverData& solver,
        QDData& QD,
        double invlambda,
        bool update_factors
){

  // compute E, C_L, C_R from the high-order solution
  if (update_factors){
    // boundary factors
    double nsumL = 0.0;
    double dsumL = 0.0;
    double nsumR = 0.0;
    double dsumR = 0.0;
    for (int u = 0; u < quadrature.order; u++){
        const double weight  = quadrature.weight[u];
        const double mu = quadrature.angle[u];
        if (mu < 0.0){
            // leaving through the left boundary
            nsumL += weight * mu * solver.psie[u].front();
            dsumL += weight * solver.psie[u].front();
        } else {
            // leaving through the right boundary
            nsumR += weight * mu * solver.psie[u].back();
            dsumR += weight * solver.psie[u].back();
        }
    }
    QD.CbfL = nsumL / dsumL;
    QD.CbfR = nsumR / dsumR;

    // Eddington factors at cell edges
    for (int i = 0; i < mesh.nedges; i++){
        double nsum = 0.0;
        double dsum = 0.0;
        for (int u = 0; u < quadrature.order; u++){
            const double mu = quadrature.angle[u];
            const double weight = quadrature.weight[u];
            nsum += weight * mu * mu * solver.psie[u][i];
            dsum += weight * solver.psie[u][i];
        }
        QD.Edde[i] = nsum / dsum;
    }

    // Eddington factors from cell averages
    for (int i = 0; i < mesh.ncells; i++){
        double nsum = 0.0;
        double dsum = 0.0;
        for (int u = 0; u < quadrature.order; u++){
            const double mu = quadrature.angle[u];
            const double weight = quadrature.weight[u];
            nsum += weight * mu * mu * solver.psibar[u][i];
            dsum += weight * solver.psibar[u][i];
        }
        QD.Eddbar[i] = nsum / dsum;
    }
  }

    // build linear system
    const int N = mesh.ncells;
    const int n = 2 * N + 1;

    std::vector<std::vector<double>> Amat(n, std::vector<double>(n, 0.0));
    std::vector<double> bvec(n, 0.0);
    std::vector<double> phivec(n, 0.0);

    // computencoefficients
    std::vector<double> a1(N);
    std::vector<double> a2(N);
    std::vector<double> b(N);
    for (int i = 0; i < N; i++){
        double sigmar = materials[mesh.mapmat[i]].xstot - materials[mesh.mapmat[i]].xsscat;
        double q = 0.0;
        if (core.search == "fixed_source"){
            sigmar -= invlambda * materials[mesh.mapmat[i]].nu * materials[mesh.mapmat[i]].xsfiss;
            q = materials[mesh.mapmat[i]].q;
        } else {
            q = qfis[i];
        }
        a1[i] = 2.0 / (materials[mesh.mapmat[i]].xstot * mesh.hx[i]);
        a2[i] = 2.0 * a1[i] * QD.Eddbar[i] + sigmar * mesh.hx[i];
        b[i] = q * mesh.hx[i];
    }

    // compute incoming flux and current
    double phiL = 0.0;
    double JinL = 0.0;
    double phiR = 0.0;
    double JinR = 0.0;
    for (int u = 0; u < quadrature.order; u++){
        const double weight = quadrature.weight[u];
        const double mu = quadrature.angle[u];
        if (mu > 0.0){
            phiL += weight * psilbc[u];
            JinL += weight * mu * psilbc[u];
        } else {
            phiR += weight * psirbc[u];
            JinR += weight * mu * psirbc[u];
        }
    }

    // left edge
    Amat[0][0] = a1[0] * QD.Edde[0];
    Amat[0][1] = -a1[0] * QD.Eddbar[0];
    if (mesh.mapbc.front() != 2){
        Amat[0][0] -= QD.CbfL;
        bvec[0] = JinL - QD.CbfL * phiL;
    }

    // balance equation
    for (int i = 0; i < N; i++){
        const int k = 2 * i + 1;
        Amat[k][k-1] = -a1[i] * QD.Edde[i];
        Amat[k][k] = a2[i];
        Amat[k][k+1] = -a1[i] * QD.Edde[i+1];
        bvec[k] = b[i];
    }

    // current continuity on interior
    for (int i = 0; i < N - 1; i++){
        const int k = 2 * i + 2;
        Amat[k][k-1] = a1[i] * QD.Eddbar[i];
        Amat[k][k] = -(a1[i] + a1[i+1]) * QD.Edde[i+1];
        Amat[k][k+1] = a1[i+1] * QD.Eddbar[i+1];
    }

    // right edge
    Amat[n-1][n-2] = a1[N-1] * QD.Eddbar[N-1];
    Amat[n-1][n-1] = -a1[N-1] * QD.Edde[N];
    if (mesh.mapbc.back() != 2){
        Amat[n-1][n-1] -= QD.CbfR;
        bvec[n-1] = JinR - QD.CbfR * phiR;
    }

    // solve linear system with the Thomas algorithm
    for (int k = 1; k < n; k++){
        const double f = Amat[k][k-1] / Amat[k-1][k-1];
        Amat[k][k] -= f * Amat[k-1][k];
        Amat[k][k-1] = 0.0;
        bvec[k] -= f * bvec[k-1];
    }
    // back substitution
    phivec[n-1] = bvec[n-1] / Amat[n-1][n-1];
    for (int k = n - 2; k >= 0; k--){
        phivec[k] = (bvec[k] - Amat[k][k+1] * phivec[k+1]) / Amat[k][k];
    }

    // update LO scalar flux values
    for (int i = 0; i < N; i++) {
        QD.phieLO[i]   = phivec[2 * i];
        QD.phibarLO[i] = phivec[2 * i + 1];
    }
    QD.phieLO[N] = phivec[2 * N];

    // update LO current values
    for (int i = 0; i < N; i++) {
        QD.jeLO[i] = a1[i] * (QD.Edde[i] * QD.phieLO[i] - QD.Eddbar[i] * QD.phibarLO[i]);
    }
    QD.jeLO[N] = a1[N-1] * (QD.Eddbar[N-1] * QD.phibarLO[N-1] - QD.Edde[N] * QD.phieLO[N]);

}