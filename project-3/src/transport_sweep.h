#ifndef TRANSPORT_SWEEP_H
#define TRANSPORT_SWEEP_H

#include "data.h"

class TransportSweep {
  public:
    void sweep_group(int g, std::vector<double> &Q, std::vector<double> &Qhat, double alpha, quadratureData &quadrature,
                     std::vector<materialData> &materials, meshData &mesh, solverData &solver);
};

#endif
