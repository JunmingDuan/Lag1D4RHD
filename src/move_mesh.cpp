#include "Lagranian1D.h"

void Lagranian1D::cal_us_roeav(Sol& ReconL_Pri, Sol& ReconR_Pri, vvector<double>& us) {
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x+1; ++i) {
    double srhol, srhor, ul, ur;
    srhol = sqrt(ReconL_Pri[i][0]);
    srhor = sqrt(ReconR_Pri[i][0]);
    ul = ReconL_Pri[i][1];
    ur = ReconR_Pri[i][1];
    us[i] = (srhol*ul + srhor*ur)/(srhol+srhor);
  }
  us[0] = 0;
  us[N_x] = 0;
}

void Lagranian1D::move_mesh(vvector<double>& mesh, vvector<double>& us, double dt, vvector<double>& mesh1) {
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x+1; ++i) {
    mesh1[i] = mesh[i] + dt * us[i];
  }
}

