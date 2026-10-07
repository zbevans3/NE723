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
        solverData& solver
    );
    
    std::vector<double> phibarold;    // old scalar flux (l-1)
    std::vector<double> phihatold;    // old first moment of scalar flux
    std::vector<double> phibarlm2;    // old scalar flux (l-2)
    std::vector<double> Qi;           // l-th iterate of fission, scattering + fixed source
    std::vector<double> Qhati;        // l-th iterate of fission, scattering first moment of source  
    std::vector<double> psilbc;       // array holding left boundary angular flux values (per ordinate)
    std::vector<double> psirbc;       // array holding right boundary angular flux values (per ordinate)
    std::vector<double> phieps;       // array of scalar flux difference

};

#endif