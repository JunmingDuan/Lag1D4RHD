#include "Lagranian1D.h"

bU Lagranian1D::LF(const bU& CONL, const bU& CONR, const bU& PRIL, const bU& PRIR, const double alpha) {
  return 0.5*(F(CONL, PRIL) + F(CONR, PRIR)) - 0.5*(CONR-CONL)*alpha;
}


void Lagranian1D::cal_flux_LF(double alpha) {
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x+1; ++i) {
    if(i == 0) FLUX[i] = LF(ghostl.Con, Con[i], ghostl.Pri, Pri[i], alpha);
    else if(i == N_x) FLUX[i] = LF(Con[i-1], ghostr.Con, Pri[i-1], ghostr.Pri, alpha);
    else FLUX[i] = LF(Con[i-1], Con[i], Pri[i-1], Pri[i], alpha);
  }
}

void Lagranian1D::forward_LF(double dt, double alpha) {
  InfiniteBD();
  cal_flux_LF(alpha);
  cal_us_roeav();
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] *= (mesh[i+1]-mesh[i]);
  }
  move_mesh(dt);
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] = (Con[i] - dt*(FLUX[i+1] - FLUX[i])) / (mesh[i+1]-mesh[i]);
  }
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Pri[i] = Con2Pri(Con[i], Gamma[i]);
  }
}


