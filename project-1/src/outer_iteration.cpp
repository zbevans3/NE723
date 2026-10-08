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

       double invlambda = solver.lambda;
       SourceIteration inner;

       std::cout << "starting outer iteration\n";
       solver.k_outer = 0;

       double rhok = 1.0e6;
       double kold = 0.0;
       double keps = 1.0e6;       
       double kepslm2 = 1.0e6;
       double sfisnew = 1.0;
       double sfisold = 0.0;
       double stopcrit = 1.0e6;
       while (solver.k_outer < core.maxout){
        solver.k_outer++;
       if (core.itdebug){
       std::cout << std::right << std::setw(12) << "Outer iteration" << std::setw(16) << solver.k_outer << "\n";
       };

        // update l-2 values
        (solver.k_outer == 1) ? kepslm2 = 1.0 : kepslm2 = keps;

        kold = solver.lambda;
        sfisold = sfisnew;

        invlambda = 1.0 / solver.lambda;

        // do inner iteration
        inner.solve(core, quadrature, materials, mesh, solver, invlambda);
         
        // compute new keff from new fission source
        sfisnew = 0.0;
        for (int i = 0 ; i < mesh.ncells ; i++){
            sfisnew += materials[mesh.mapmat[i]].nu * materials[mesh.mapmat[i]].xsfiss * solver.phibar[i];
        }
        solver.lambda = kold * (sfisnew / sfisold);

        keps = fabs(solver.lambda-kold);
        rhok = keps / kepslm2;

        if (core.itdebug){
            std::cout << std::right << "keff " << std::setw(16) << std::fixed << std::setprecision(8) << solver.lambda
                      << std::setw(8) << "error" << std::setw(16) << std::scientific << std::setprecision(6) << keps << '\n';
        }

        // compute rhok
        stopcrit = core.epsk * (1.0 / rhok - 1.0);
        if (keps < stopcrit) {
            std::cout << "Outer iteration convergence reached in " << solver.k_outer << " iterations.\n";
            std::cout << "Final convergence criterion: " 
                      << std::scientific << std::setprecision(6) << stopcrit << "\n";
            std::cout << "\n";
            break;
        }
        std::cout << "\n";

       }

}