#ifndef SOURCE_ITERATION_H
#define SOURCE_ITERATION_H

#include "data.h"

class SourceIteration {
  public:
    void solve(coreData &core, quadratureData &quadrature, std::vector<materialData> &materials, meshData &mesh,
               solverData &solver, QDData &QD, QDMGData &QDMG, QD1GData &QD1G, double invlambda);

  private:
    void compute_residuals(coreData &core, std::vector<materialData> &materials, meshData &mesh, solverData &solver,
                           double &alpha, double invlambda);

    void update_flux_and_current(quadratureData &quadrature, meshData &mesh, solverData &solver, int g);

    void QD_acceleration(coreData &core, quadratureData &quadrature, std::vector<materialData> &materials,
                         meshData &mesh, solverData &solver, QDData &QD, double invlambda, bool update_factors);

    std::vector<std::vector<double>> phibarold; // old scalar flux (l-1)
    std::vector<std::vector<double>> phihatold; // old first moment of scalar flux
    std::vector<std::vector<double>> Qi;        // l-th iterate of fission, scattering + fixed source
    std::vector<std::vector<double>> Qhati;     // l-th iterate of fission, scattering first moment of source
    std::vector<std::vector<double>> phieps;    // array of scalar flux difference
    std::vector<std::vector<double>> qfis;      // array of fission source
    std::vector<std::vector<double>> qfishat;   // array of first moment fission source
    std::vector<double> JLOe;                   // low-order edge current (nedges)
    std::vector<std::vector<double>> qupscat;   // upscattering source
    std::vector<double> Al;
    std::vector<double> Ad;
    std::vector<double> Au;
    std::vector<double> bvec;
    std::vector<double> phivec;
    std::vector<double> a1;
    std::vector<double> a2;
    std::vector<double> b;
};

#endif