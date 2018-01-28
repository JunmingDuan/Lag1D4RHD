#include "Lagranian1D.h"

void Lagranian1D::cal_us_roeav(Sol& ReconL_Pri, Sol& ReconR_Pri,
    const VEC& GAMMAL, const VEC& GAMMAR, vvector<double>& us) {
  // gamma*u = <sqrt(\rho h)*gamma*u> / <sqrt(\rho h)>, <.> means arithmetic mean
  // gamma = <sqrt(\rho h)*gamma> / <sqrt(\rho h)>
  // u = gamma*u / gamma;
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x+1; ++i) {
    double hl, hr, kl, kr, gl, gr, ul, ur;
    hl = 1 + GAMMAL[i]/(GAMMAL[i]-1)*ReconL_Pri[i][2]/ReconL_Pri[i][0];
    hr = 1 + GAMMAR[i]/(GAMMAR[i]-1)*ReconR_Pri[i][2]/ReconR_Pri[i][0];
    kl = sqrt(ReconL_Pri[i][0]*hl);
    kr = sqrt(ReconR_Pri[i][0]*hr);
    ul = ReconL_Pri[i][1];
    ur = ReconR_Pri[i][1];
    gl = 1./sqrt(1-ul*ul);
    gr = 1./sqrt(1-ur*ur);
    double tmp1 = (kl*gl*ul + kr*gr*ur)/(kl+kr);
    double tmp2 = (kl*gl + kr*gr)/(kl+kr);
    us[i] = tmp1/tmp2;
    //us[i] = (sqrt(ReconL_Pri[i][0])*ul + sqrt(ReconR_Pri[i][0])*ur) / (sqrt(ReconL_Pri[i][0]) + sqrt(ReconR_Pri[i][0]));
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

