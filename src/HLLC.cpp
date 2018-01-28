#include "Lagranian1D.h"

bU Lagranian1D::HLLC(const bU& CONL, const bU& CONR, const bU& PRIL, const bU& PRIR, const double Gammal, const double Gammar, double& SM) {
  double SL, SR;//, ui, ci, vl, vr;
  double hl, csl, ul, hr, csr, ur;
  double lam1, lam2, lam3, lam4;

  //characteristic speed of left and right side
  hl = 1+PRIL[2]/PRIL[0]*Gammal/(Gammal-1);
  csl = sqrt(Gammal*PRIL[2]/PRIL[0]/hl);
  ul = PRIL[1];
  lam1 = ul-(1-ul*ul)*csl/(1-ul*csl);//negative speed
  lam2 = ul+(1-ul*ul)*csl/(1+ul*csl);//positive speed
  hr = 1+PRIR[2]/PRIR[0]*Gammar/(Gammar-1);
  csr = sqrt(Gammar*PRIR[2]/PRIR[0]/hr);
  ur = PRIR[1];
  lam3 = ur-(1-ur*ur)*csr/(1-ur*csr);//negative speed
  lam4 = ur+(1-ur*ur)*csr/(1+ur*csr);//positive speed
  SL = std::min((lam1), (lam3));//left characteristic speed
  SR = std::max((lam2), (lam4));//right characteristic speed
  {
    //cal S^* and p^*, S^* is the smaller root of a quadric equation:
    //F_E^{hll} x^2 - (E^{hll}+F^{hll}_m) x + m^{hll} = 0
    double coe1 = SR*CONL[1] - SL*CONR[1] + SL*SR*(CONR[2] - CONL[2]);
    double coe2 = ( SR*(CONL[1]*PRIL[1]+PRIL[2]) - SL*(CONR[1]*PRIR[1]+PRIR[2]) + SL*SR*(CONR[1] - CONL[1])
        + SR*CONR[2] - SL*CONL[2] + CONL[1] - CONR[1] );
    double coe3 = SR*CONR[1] - SL*CONL[1] + CONL[1]*PRIL[1]+PRIL[2] - CONR[1]*PRIR[1]-PRIR[2];
    double PM;
    if(fabs(coe1) < 1e-15) SM = coe3/coe2;
    else if(fabs(coe1*coe3) < 1e-15) SM = 0;
    else SM = (coe2 - sqrt(coe2*coe2 - 4.*coe1*coe3))/2./coe1;
    PM = (SM*(SL*CONL[2]-CONL[1]) + PRIL[2] - CONL[1]*(SL-PRIL[1])) / (1-SL*SM);
    //std::cout.setf(std::ios::scientific);
    //std::cout.precision(16);
    //std::cout << coe1 << " " << coe2 << " " << coe3 << " " << (SM*(SL*CONL[2]-CONL[1]) + PRIL[2] - CONL[1]*(SL-PRIL[1])) / (1-SL*SM) -
    //(SM*(SR*CONR[2]-CONR[1]) + PRIR[2] - CONR[1]*(SR-PRIR[1])) / (1-SR*SM) << std::endl;
    bU F;
    F[0] = 0;
    F[1] = PM;
    F[2] = PM*SM;
    return F;
  }
}

void Lagranian1D::cal_flux_HLLC(Sol& ReconL_Con, Sol& ReconR_Con, Sol& ReconL_Pri, Sol& ReconR_Pri,
    Sol& FLUX, vvector<double>& us) {
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x+1; ++i) {
    if(i == 0) FLUX[i] = HLLC(ReconL_Con[i], ReconR_Con[i], ReconL_Pri[i], ReconR_Pri[i],
        Gamma[0], Gamma[i], us[i]);
    else if(i == N_x) FLUX[i] = HLLC(ReconL_Con[i], ReconR_Con[i], ReconL_Pri[i], ReconR_Pri[i],
        Gamma[i-1], Gamma[N_x-1], us[i]);
    else FLUX[i] = HLLC(ReconL_Con[i], ReconR_Con[i], ReconL_Pri[i], ReconR_Pri[i],
        Gamma[i-1], Gamma[i], us[i]);
  }
}

