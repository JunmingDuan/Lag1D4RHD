#include "Lagranian1D.h"

bU Lagranian1D::LLF(const bU& CONL, const bU& CONR, const bU& PRIL, const bU& PRIR, const double alpha) {
  return 0.5*(F(CONL, PRIL) + F(CONR, PRIR)) - 0.5*(CONR-CONL)*alpha;
}

void Lagranian1D::cal_flux_LLF(Sol& ReconL_Con, Sol& ReconR_Con, Sol& ReconL_Pri, Sol& ReconR_Pri, Sol& FLUX) {
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x+1; ++i) {
    double alpha, laml, lamr, roe_lam1, roe_lam3;
    double roe_lam;
    //cal local characteristic speed
    if(i == 0) {
      laml = cal_max_lambda_Lag(ReconL_Con[i], ReconL_Pri[i], Gamma[0]);
      lamr = cal_max_lambda_Lag(ReconR_Con[i], ReconR_Pri[i], Gamma[i]);
      //cal_min_max_roe_lam(ReconL_Con[i], ReconR_Con[i], ReconL_Pri[i], ReconR_Pri[i], Gamma[0], Gamma[i], roe_lam1, roe_lam3);
      roe_lam = cal_max_roe_lam_Lag(ReconL_Con[i], ReconR_Con[i], ReconL_Pri[i], ReconR_Pri[i], Gamma[0], Gamma[i]);
    }
    else if(i == N_x) {
      laml = cal_max_lambda_Lag(ReconL_Con[i], ReconL_Pri[i], Gamma[i-1]);
      lamr = cal_max_lambda_Lag(ReconR_Con[i], ReconR_Pri[i], Gamma[N_x-1]);
      //cal_min_max_roe_lam(ReconL_Con[i], ReconR_Con[i], ReconL_Pri[i], ReconR_Pri[i], Gamma[i-1], Gamma[N_x-1], roe_lam1, roe_lam3);
      roe_lam = cal_max_roe_lam_Lag(ReconL_Con[i], ReconR_Con[i], ReconL_Pri[i], ReconR_Pri[i], Gamma[i-1], Gamma[N_x-1]);
    }
    else {
      laml = cal_max_lambda_Lag(ReconL_Con[i], ReconL_Pri[i], Gamma[i-1]);
      lamr = cal_max_lambda_Lag(ReconR_Con[i], ReconR_Pri[i], Gamma[i]);
      //cal_min_max_roe_lam(ReconL_Con[i], ReconR_Con[i], ReconL_Pri[i], ReconR_Pri[i], Gamma[i-1], Gamma[i], roe_lam1, roe_lam3);
      roe_lam = cal_max_roe_lam_Lag(ReconL_Con[i], ReconR_Con[i], ReconL_Pri[i], ReconR_Pri[i], Gamma[i-1], Gamma[i]);
    }
    //alpha = std::max(std::max(laml, lamr), std::max(fabs(roe_lam1), fabs(roe_lam3)));
    alpha = std::max(std::max(laml, lamr), roe_lam);
    //alpha = std::max(laml, lamr);
    FLUX[i] = LLF(ReconL_Con[i], ReconR_Con[i], ReconL_Pri[i], ReconR_Pri[i], 1.0*alpha);
  }
}

