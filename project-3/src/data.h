#ifndef DATA_H
#define DATA_H

#include <string>
#include <vector>

struct coreData {
    int nregions = 0;          // number of regions
    std::vector<int> nx;       // cells per region
    std::vector<double> xedge; // right edge of cell (most left edge is assumed to be 0.0)
    int nmat = 0;              // number of unique materials
    std::vector<int> matid;    // material ID of material in region
    int ngroups = 1;           // number of energy groups

    int maxin = 100;        // maximum number of inner (SI) iterations
    int maxout = 50;        // maximum number of outer (eigenvalue) iterations
    double epsk = 1.0e-6;   // convergence criterion on eigenvalue
    double epsflx = 1.0e-6; // convergence criterion on flux
    double epspow = 1.0e-6; // convergence criterion on fission source
    double kcrit = 1.0;     // critical eigenvalue for critical size search

    std::string search = "fixed_source"; // solution type; (fixed_source, eigenvalue, critical_size)
    std::string scheme = "LD";           // discretization scheme; (LD, LLD)
    std::string bcleft = "vacuum";       // left boundary condition; (vacuum, reflective, incoming)
    std::string bcright = "vacuum";      // right boundary condition; (vacuum, reflective, incoming)
    std::string acceleration = "none";   // acceleration scheme; (none, QD, MLQD)
    double psileft = 0.0;                // angular flux at left boundary
    double psiright = 0.0;               // angular flux at right boundary

    bool itdebug = true; // option to print inner iteration data
    bool otdebug = true; // option to print outer iteration data
};

struct meshData {
    int ncells = 0;              // total number of cells
    int nedges = 0;              // total number of edges (ncells + 1)
    std::vector<int> mapmat;     // material ID in cell i mapmat(ncells)
    std::vector<int> mapbc;      // boundary condition in at edge i mapbc(nedges)
    std::vector<double> hx;      // i-th cell width hx(ncells)
    std::vector<double> xloc;    // x-location of i-th edge x(nedges)
    std::vector<double> xcenter; // x-location of i-th center x(ncells)
};

struct solverData {
    std::vector<std::vector<std::vector<double>>> psie;   // cell edge angular flux (ngroups, nangles, nedges)
    std::vector<std::vector<std::vector<double>>> psibar; // cell average angular flux (ngroups, nangles, ncells)
    std::vector<std::vector<std::vector<double>>> psihat; // first moment of angular flux (ngroups, nangles, ncells)
    std::vector<std::vector<double>> psilbc;              // left BC edge angular flux (ngroups, nangles)
    std::vector<std::vector<double>> psirbc;              // right BC edge angular flux (ngroups, nangles)
    std::vector<std::vector<double>> phie;                // cell edge scalar flux (ngroups, nedges)
    std::vector<std::vector<double>> phibar;              // cell average scalar flux (ngroups, ncells)
    std::vector<std::vector<double>> phihat;              // first moment of scalar flux (ngroups, ncells)
    std::vector<std::vector<double>> je;                  // cell edge current (ngroups, nedges)
    std::vector<std::vector<double>> jbar;                // cell average current (ngroups, ncells)
    std::vector<std::vector<double>> jhat;                // first moment of current (ngroups, ncells)
    std::vector<std::vector<double>> source;              // RHS of zeroth moment eqm (ngroups, ncells)
    std::vector<std::vector<double>> source1;             // RHS of first moment eqn (ngroups, ncells)
    std::vector<std::vector<double>> sresidual;           // scheme residual for zeroth moment eqn (ngroups, ncells)
    std::vector<std::vector<double>> sresidual1;          // scheme residual for first moment eqn (ngroups, ncells)
    std::vector<std::vector<double>> iresidual;           // iterative residual for zeroth moment eqn (ngroups, ncells)
    std::vector<std::vector<double>> iresidual1;          // iterative residual for first moment eqn (ngroups, ncells)
    std::vector<double> ispecrad;                         // specral radius vs inner iterate
    std::vector<double> ilinfphi;                         // infinity norm of scalar flux vs inner iterate
    std::vector<double> ospecrad;                         // specral radius vs inner iterate
    std::vector<double> olinfphi;                         // infinity norm of scalar flux vs inner iterate
    std::vector<double> specradk;                         // specral radius vs iterate of keff
    std::vector<double> linfkeff;                         // infinity norm of keff flux vs iterate
    std::vector<double> keff;                             // keff flux vs iterate
    std::vector<double> linffsrc;                         // infinity norm of fissio source vs iterate
    int k_inner = 0;
    int k_outer = 0;
    double lambda = 1.0;
    double fluxnorm = 0.0;
};

struct QDData {
    std::vector<double> phieLO;   // cell edge scalar flux from low order solution
    std::vector<double> phibarLO; // cell average scalar flux from low order solution
    std::vector<double> jeLO;     // cell edge current from low order solution
    double CbfL;                  // Left boundary factor
    double CbfR;                  // Right boundary factor
    std::vector<double> Edde;     // Cell edge Eddington factor
    std::vector<double> Eddbar;   // Cell average Eddington factor
    std::vector<double> mLO;      // cell average value of the metric of transport effects
    std::vector<double> mprimeLO; // cell average value of the derivative based metric of transport effects
};

struct QDMGData {
    std::vector<std::vector<double>> Eddbar; // Cell average group Eddington factor
    std::vector<std::vector<double>> Edde;   // Cell edge group Eddington factor
    std::vector<double> CbfL;                // Left boundary group factor
    std::vector<double> CbfR;                // Right boundary group factor
    std::vector<std::vector<double>> phibar; // cell average group scalar flux from low order solution
    std::vector<std::vector<double>> phie;   // cell edge group scalar flux from low order solution
    std::vector<std::vector<double>> je;     // cell edge group current from low order solution
};

struct QD1GData {
    std::vector<double> phibar; // cell average scalar flux from one-group low order solution, including boundary values
    std::vector<double> phie;   // cell edge scalar flux from one-group low order solution
    std::vector<double> je;     // cell edge current from one-group low order solution
    std::vector<double> xsabs;  // cell average one-group absorption cross section
    std::vector<double> xsnusigf; // cell average one-group nu times fission cross section
    std::vector<double> Hm;       // cell edge H- factor
    std::vector<double> Hp;       // cell edge H+ factor
    std::vector<double> Am;       // cell edge A- factor
    std::vector<double> Ap;       // cell edge A+ factor
    double CbfL;                  // Left boundary factor
    double CbfR;                  // Right boundary factor
    double k = 1.0;               // one-group eigenvalue
};

struct materialData {
    int id = 0;
    std::vector<double> xstot;
    std::vector<double> xsabs;
    std::vector<double> xsscat;                 // within-group scattering (g -> g)
    std::vector<std::vector<double>> xsscatmat; // full scattering matrix [g][gp] = g -> gp
    std::vector<double> xsfiss;
    std::vector<double> xschi;
    std::vector<double> nu;
    std::vector<double> q;
};

struct quadratureData {
    int order = 0;
    std::vector<double> angle;
    std::vector<double> weight;
    std::vector<int> reflmap;
};

#endif
