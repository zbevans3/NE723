#ifndef DATA_H
#define DATA_H

#include <string>
#include <vector>

struct coreData {
    int nregions = 0;                // number of regions
    std::vector<int> nx;             // cells per region
    std::vector<double> xedge;       // right edge of cell (most left edge is assumed to be 0.0)
    int nmat = 0;                    // number of unique materials
    std::vector<int> matid;          // material ID of material in region

    int maxin = 100;                 // maximum number of inner (SI) iterations
    int maxout = 50;                 // maximum number of outer (eigenvalue) iterations 
    double epsk = 1.0e-6;            // convergence criterion on eigenvalue
    double epsflx = 1.0e-6;          // convergence criterion on flux

    std::string search;
    std::string bcleft = "vacuum";   // left boundary condition; (vacuum, reflective, incoming)
    std::string bcright = "vacuum";  // right boundary condition; (vacuum, reflective, incoming)
    double psileft = 0.0;            // angular flux at left boundary
    double psiright = 0.0;           // angular flux at right boundary
};

struct meshData {
    int ncells = 0;               // total number of cells
    int nedges = 0;               // total number of edges (ncells + 1)
    std::vector<int> mapmat;      // material ID in cell i mapmat(ncells)
    std::vector<int> mapbc;       // boundary condition in at edge i mapbc(nedges)
    std::vector<double> hx;       // i-th cell width hx(ncells)
    std::vector<double> xloc;        // x-location of i-th edge x(nedges)
};

struct solverData {
    std::vector<std::vector<double>> psie;
    std::vector<std::vector<double>> psim;
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
    int order = 0;
    std::vector<double> angle;
    std::vector<double> weight;
};

#endif
