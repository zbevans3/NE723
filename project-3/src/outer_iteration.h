#ifndef OUTER_ITERATION_H
#define OUTER_ITERATION_H

#include "data.h"

class OuterIteration {
  public:
    void solve(coreData &core, quadratureData &quadrature, std::vector<materialData> &materials, meshData &mesh,
               solverData &solver, QDData &QD, QDMGData &QDMG, QD1GData &QD1G);

    std::vector<double> Al;
    std::vector<double> Ad;
    std::vector<double> Au;
    std::vector<double> bvec;
    std::vector<double> phivec;
    std::vector<double> a1;
    std::vector<double> a2;
    std::vector<double> a3;
    std::vector<double> b;
    std::vector<double> hxe;
    std::vector<double> phi1gold;
};

#endif