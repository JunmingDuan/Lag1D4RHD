#include "Lagranian1D.h"

void Lagranian1D::Reconstruction(const Sol& sol, const bU& gl, const bU& gr,
    Sol& ReconL_Con, Sol& ReconR_Con, Sol& ReconL_Pri, Sol& ReconR_Pri) {
  if(is_RECON == 0) {
#pragma omp parallel for num_threads(Nthread)
    for(int i = 0; i < N_x; ++i) {
      ReconL_Con[i+1] = sol[i];
      ReconR_Con[i] = sol[i];
      ReconL_Pri[i+1] = Con2Pri(ReconL_Con[i+1], Gamma[i]);
      ReconR_Pri[i] = Con2Pri(ReconR_Con[i], Gamma[i]);
    }
    ReconL_Con[0] = gl;
    ReconR_Con[N_x] = gr;
    ReconL_Pri[0] = Con2Pri(ReconL_Con[0]);
    ReconR_Pri[N_x] = Con2Pri(ReconR_Con[N_x]);
  }
}

