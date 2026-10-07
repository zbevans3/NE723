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
            mesh.mapmat[ctr] = core.matid[ii]-1;
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
      3 - incoming
    */
    if (core.bcleft == "vacuum"){
        mesh.mapbc.front() = 1;
    } else if (core.bcleft == "reflective"){
        mesh.mapbc.front() = 2;
    } else {
        mesh.mapbc.front() = 3;
    }
     
    if (core.bcright == "vacuum"){
        mesh.mapbc.back() = 1;
    } else if (core.bcright == "reflective"){
        mesh.mapbc.back() = 2;
    } else {
        mesh.mapbc.back() = 3;
    }
}