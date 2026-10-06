#include "build_mesh.h"

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <string>

void Mesh::build(
    coreData& core,
    meshData& mesh)
{

    mesh.ncells = std::accumulate(core.nx.begin(), core.nx.end(), 0);
    mesh.nedges = mesh.ncells + 1;

    mesh.mapmat.assign(mesh.ncells, 0);
    mesh.hx.assign(mesh.ncells, 0.0);
    mesh.mapbc.assign(mesh.nedges, 0);    // defualt all nodes to interior
    mesh.xloc.assign(mesh.nedges, 0.0);

    int ctr = 0;
    double xleft = 0.0;

    // build mapmat, cell width, and cell edge x-locations
    for (int ii = 0; ii < core.nregions; ++ii) {

        double xright = core.xedge[ii];
        double dx = (xright - xleft) / static_cast<double>(core.nx[ii]);

        for (int j = 0; j < core.nx[ii]; ++j) {
            mesh.mapmat[ctr] = core.matid[ii];
            mesh.hx[ctr] = dx;
            mesh.xloc[ctr + 1] = mesh.xloc[ctr] + dx;
            ++ctr;
        }
        xleft = xright;
    }

    /* assign BC to edges
      0 - interior
      1 - vacuum
      2 - reflective
    */
    mesh.mapbc.front() = (core.bcleft == "vacuum")  ? 1 : 2; 
    mesh.mapbc.back()  = (core.bcright == "vacuum") ? 1 : 2;
}