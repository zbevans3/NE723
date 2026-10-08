#include "outer_iteration.h"

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <string>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <iomanip>

#include "source_iteration.h"

void OuterIteration::solve(
        coreData& core,
        quadratureData& quadrature,
        std::vector<materialData>& materials,
        meshData& mesh,
        solverData& solver
        )    
{

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
       std::vector<double> phiold;
       std::vector<double> sfisnew;
       std::vector<double> sfisold;
       std::vector<double> epsflx;
       std::vector<double> sfiseps;

       phiold.assign(mesh.ncells, 1.0);
       epsflx.assign(mesh.ncells, 0.0);
       sfisnew.assign(mesh.ncells, 1.0);
       sfisold.assign(mesh.ncells, 1.0);       
       sfiseps.assign(mesh.ncells, 0.0);       
       
       double invlambda = solver.lambda;
       SourceIteration inner;

       if (core.search != "critical_size") std::cout << "starting outer iteration\n";
       solver.k_outer = 0;

        if (core.otdebug){
        std::cout << std::right
             << std::setw(20) << "Outer iteration"
             << std::setw(16) << "keff"
             << std::setw(16) << "keff error"
             << std::setw(16) << "flux error"
             << std::setw(16) << "fsource error"
             << '\n';
        }

       int lout = 0;
       while (lout < core.maxout){
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
            inner.solve(core, quadrature, materials, mesh, solver, invlambda);          
    
            // compute new keff from new fission source
            sfistotnew = 0.0;
            for (int i = 0 ; i < mesh.ncells ; i++){
                sfisnew[i] = materials[mesh.mapmat[i]].nu * materials[mesh.mapmat[i]].xsfiss * solver.phibar[i];
                sfistotnew += sfisnew[i] * mesh.hx[i];
                sfiseps[i] = std::fabs(sfisnew[i] - sfisold[i]) / std::max(std::fabs(sfisold[i]), 1.0e-14);
                epsflx[i] = std::fabs(solver.phibar[i] - phiold[i]);
            }
            maxsfeps = *std::max_element(sfiseps.begin(), sfiseps.end());
            maxepsflx = *std::max_element(epsflx.begin(), epsflx.end());
            rhoflx = maxepsflx / maxepsflxlm2;

            solver.lambda = kold * (sfistotnew / sfistotold);
            solver.keff.push_back(solver.lambda);

            keps = fabs(solver.lambda-kold);
            rhok = keps / kepslm2;

            solver.linfkeff.push_back(keps);
            solver.specradk.push_back(rhok);
            solver.olinfphi.push_back(maxepsflx);
            solver.ospecrad.push_back(rhoflx);
            solver.linffsrc.push_back(maxsfeps);

             if (core.otdebug){
             std::cout << std::right
                  << std::setw(20) << lout
                  << std::setw(16) << std::fixed << std::setprecision(8) << solver.lambda
                  << std::setw(16) << std::scientific << std::setprecision(6) << keps
                  << std::setw(16) << std::scientific << std::setprecision(6) << maxepsflx                  
                  << std::setw(16) << std::scientific << std::setprecision(6) << maxsfeps
                  << '\n';
             }
    
            // compute stopping criterion
            stopkcrit = core.epsk * (1.0 / rhok - 1.0);
            stopflx = core.epsflx * (1.0 / rhoflx - 1.0);
            if (keps < stopkcrit && stopflx < core.epsflx && maxsfeps < core.epspow) {
                if (core.search != "critical_size"){
                   std::cout << "Outer iteration convergence reached in " << lout << " iterations.\n";
                   std::cout << "Final convergence criterion: \n";
                   std::cout << "  keff" << std::setw(16) << std::scientific << std::setprecision(6) << stopkcrit << "\n";
                   std::cout << "  flux" << std::setw(16) << std::scientific << std::setprecision(6) << stopflx << "\n";                   
                   std::cout << "\n";
                }
                break;            
            }
       }

       // normalize flux to 1.0
       double phitot = std::accumulate(solver.phibar.begin(), solver.phibar.end(), 0.0);
       for (int i = 0 ; i < mesh.ncells ; i++){
           solver.phibar[i] /= phitot;
       }

       solver.k_outer += lout;

}