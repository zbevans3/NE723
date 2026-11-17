#include "transport_sweep.h"

void TransportSweep::sweep_group(int g,
                 std::vector<double>& Q,
                 std::vector<double>& Qhat,
                 double alpha,
                 quadratureData& quadrature,
                 std::vector<materialData>& materials,
                 meshData& mesh,
                 solverData& solver)
{
    std::vector<std::vector<double>>& psie = solver.psie[g];
    std::vector<std::vector<double>>& psibar = solver.psibar[g];
    std::vector<std::vector<double>>& psihat = solver.psihat[g];

    for (int u = 0; u < quadrature.order; u++){
        const double mu = quadrature.angle[u];
        const double invmu = 1.0 / mu;
        if (mu > 0.0){
            psie[u].front() = solver.psilbc[g][u];
            for (int i = 0; i < mesh.ncells; i++){
                const double betai = materials[mesh.mapmat[i]].xstot[g] * mesh.hx[i] * invmu;
                const double invGi = 1.0 / (alpha*betai*betai + (3.0+alpha)*betai + 6.0);
                psibar[u][i] = ((3.0 + alpha*betai)*Q[i] - alpha*Qhat[i]) * invmu * invGi
                             + (6.0 + alpha*betai) * invGi * psie[u][i];
                psihat[u][i] = (3.0*Q[i] + alpha*(1.0 + betai)*Qhat[i]) * invmu * invGi
                             - 3.0*betai * invGi * psie[u][i];
                psie[u][i+1] = psibar[u][i] + psihat[u][i];
            }
        } else {
            psie[u].back() = solver.psirbc[g][u];
            for (int i = mesh.ncells - 1; i >= 0; --i){
                const double betai = materials[mesh.mapmat[i]].xstot[g] * mesh.hx[i] * invmu;
                const double invGi = 1.0 / (alpha*betai*betai - (3.0+alpha)*betai + 6.0);
                psibar[u][i] = ((alpha*betai - 3.0)*Q[i] - alpha*Qhat[i]) * invmu * invGi
                             + (6.0 - alpha*betai) * invGi * psie[u][i+1];
                psihat[u][i] = (3.0*Q[i] + alpha*(betai - 1.0)*Qhat[i]) * invmu * invGi
                             - 3.0*betai * invGi * psie[u][i+1];
                psie[u][i] = psibar[u][i] - psihat[u][i];
            }
        }
    }

    for (int u = 0; u < quadrature.order; ++u) {
        if (mesh.mapbc.front() == 2 && quadrature.angle[u] > 0.0) {
            solver.psilbc[g][u] = psie[quadrature.reflmap[u]].front();
        }
        if (mesh.mapbc.back() == 2 && quadrature.angle[u] < 0.0) {
            solver.psirbc[g][u] = psie[quadrature.reflmap[u]].back();
        }
    }

}