#ifndef DATA_H
#define DATA_H

#include <string>
#include <vector>

struct coreData {
    int nsize = 0;
    int isize = 0;
    int maxin = 0;
    int maxout = 0;
    double hx = 0.0;
    double epsk = 0.0;
    double epspow = 0.0;
    int nmat = 0;
    std::string search;
    std::vector<int> mapbc;
    std::vector<int> mapmat;
};

struct solverData {
    std::vector<std::vector<double>> psi;
    std::vector<std::vector<double>> source;
    std::vector<double> residual;
    int k_inner = 0;
    int k_outer = 0;
    bool kFlag = false;
    double lambda = 1.0;
    double fluxnorm = 0.0;
    int iflag = 0;
};

struct materialData {
    int id = 0;
    double xstot = 0.0;
    double xsscat = 0.0;
    double xsfiss = 0.0;
    double nu = 0.0;
    double q = 0.0;
};

struct quadratureData {
    int nordinates = 0;
    std::vector<double> angle;
    std::vector<double> weight;
};

#endif
