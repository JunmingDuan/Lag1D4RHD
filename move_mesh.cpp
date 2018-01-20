#include "Lagranian1D.h"

void Lagranian1D::cal_us_roeav() {
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x+1; ++i) {
    double srhol, srhor, ul, ur;
    if(i == 0) {
      srhol = sqrt(ghostl.Pri[0]);
      ul = sqrt(ghostl.Pri[1]);
    }
    else {
      srhol = sqrt(Pri[i-1][0]);
      ul = Pri[i-1][1];
    }
    if(i == N_x) {
      srhor = sqrt(ghostr.Pri[0]);
      ur = sqrt(ghostr.Pri[1]);
    }
    else {
      srhor = sqrt(Pri[i][0]);
      ur = Pri[i][1];
    }
    us[i] = (srhol*ul + srhor*ur)/(srhol+srhor);
  }
  us[0] = 0;
  us[N_x] = 0;
}

void Lagranian1D::move_mesh(double dt) {
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x+1; ++i) {
    mesh[i] += dt * us[i];
  }
}

