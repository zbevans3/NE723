#ifndef SOURCE_ITERATION_H
#define SOURCE_ITERATION_H

#include "data.h"

class SourceIteration {
public:
    void solve(
        coreData& core,
        quadratureData& quadrature,
        std::vector<materialData>& materials,
        meshData& mesh,
        solverData& solver,
        double invlambda
    );
    
private:
    void update_flux_and_current(
            quadratureData& quadrature,
            meshData& mesh,
            solverData& solver
    );

    void QD_acceleration(
            coreData& core,
            quadratureData& quadrature,
            std::vector<materialData>& materials,
            meshData& mesh,
            solverData& solver,
            double invlambda
    );

    std::vector<double> phibarold;    // old scalar flux (l-1)
    std::vector<double> phihatold;    // old first moment of scalar flux
    std::vector<double> Qi;           // l-th iterate of fission, scattering + fixed source
    std::vector<double> Qhati;        // l-th iterate of fission, scattering first moment of source  
    std::vector<double> psilbc;       // array holding left boundary angular flux values (per ordinate)
    std::vector<double> psirbc;       // array holding right boundary angular flux values (per ordinate)
    std::vector<double> phieps;       // array of scalar flux difference
    std::vector<double> qfis;         // array of fission source
    std::vector<double> qfishat;      // array of first moment fission source

    // QD data
    double CbfL;         // left edge boundary factor
    double CbfR;         // right edge boundary factor
    std::vector<double> Edde;         // Edge Eddington factor
    std::vector<double> Eddbar;       // Cell average Eddington factor
    std::vector<double> phiLOe;       // low-order edge scalar flux (nedges)
    std::vector<double> phiLObar;     // low-order cell scalar flux (ncells)
    std::vector<double> JLOe;         // low-order edge current (nedges)
};

#endif