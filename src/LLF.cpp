#include "Lagranian1D.h"

bU Lagranian1D::LLF(const bU& CONL, const bU& CONR, const bU& PRIL, const bU& PRIR, const double alpha) {
  return 0.5*(F(CONL, PRIL) + F(CONR, PRIR)) - 0.5*(CONR-CONL)*alpha;
}

void Lagranian1D::cal_flux_LLF(Sol& Con, Sol& Pri, Sol& FLUX) {
  double alpha;
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x+1; ++i) {
    if(i == 0) {
      //cal local characteristic speed
      double laml = cal_max_lambda_Lag(ghostl.Con, ghostl.Pri, ghostl.Gamma);
      double lamr = cal_max_lambda_Lag(i);
      alpha = std::max(laml, lamr);
      FLUX[i] = LLF(ghostl.Con, Con[i], ghostl.Pri, Pri[i], alpha);
    }
    else if(i == N_x) {
      double laml = cal_max_lambda_Lag(i-1);
      double lamr = cal_max_lambda_Lag(ghostr.Con, ghostr.Pri, ghostr.Gamma);
      alpha = std::max(laml, lamr);
      FLUX[i] = LLF(Con[i-1], ghostr.Con, Pri[i-1], ghostr.Pri, alpha);
    }
    else {
      double laml = cal_max_lambda_Lag(i-1);
      double lamr = cal_max_lambda_Lag(i);
      alpha = std::max(laml, lamr);
      FLUX[i] = LLF(Con[i-1], Con[i], Pri[i-1], Pri[i], alpha);
    }
  }
}

