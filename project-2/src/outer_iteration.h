#ifndef OUTER_ITERATION_H
#define OUTER_ITERATION_H

#include "data.h"

class OuterIteration {
public:
    void solve(
        coreData& core,
        quadratureData& quadrature,
        std::vector<materialData>& materials,
        meshData& mesh,
        solverData& solver
    );

};

#endif