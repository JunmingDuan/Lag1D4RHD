#include "Lagranian1D.h"

bU Lagranian1D::HLLC(bU& CONL, bU& CONR, bU& PRIL, bU& PRIR, const double Gammal, const double Gammar, double& SM) {
  double SL, SR;//, ui, ci, vl, vr;
  double hl, csl, ul, hr, csr, ur;
  double lam1, lam2, lam3, lam4;
  double roe_lam1, roe_lam3;

  //characteristic speed of left and right side
  cal_min_max_roe_lam(CONL, CONR, PRIL, PRIR, Gammal, Gammar, roe_lam1, roe_lam3);
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
  //SL = std::min(lam1, lam3);//left characteristic speed
  //SR = std::max(lam2, lam4);//right characteristic speed
  SL = std::min(std::min(lam1, lam3), roe_lam1);//left characteristic speed
  SR = std::max(std::max(lam2, lam4), roe_lam3);//right characteristic speed
  //bU CHARL, CHARR;
  //std::cout << "before\n" << CONL << "\n" << CONR << std::endl;
  //CHAR_DECOM(CONL, CONR, PRIL, PRIR, Gammal, Gammar, CHARL, CHARR, 1);
  //std::cout << "after\n" << CHARL << "\n" << CHARR << std::endl;
  //CHAR_DECOM(CONL, CONR, PRIL, PRIR, Gammal, Gammar, CHARL, CHARR, -1);
  //std::cout << "back\n" << CONL << "\n" << CONR << std::endl;
  {
    //cal S^* and p^*, S^* is the smaller root of a quadric equation:
    //F_E^{hll} x^2 - (E^{hll}+F^{hll}_m) x + m^{hll} = 0
    double coe1 = SR*CONL[1] - SL*CONR[1] + SL*SR*(CONR[2] - CONL[2]);
    double coe2 = ( SR*(CONL[1]*PRIL[1]+PRIL[2]) - SL*(CONR[1]*PRIR[1]+PRIR[2]) + SL*SR*(CONR[1] - CONL[1])
        + SR*CONR[2] - SL*CONL[2] + CONL[1] - CONR[1] );
    double coe3 = SR*CONR[1] - SL*CONL[1] + CONL[1]*PRIL[1]+PRIL[2] - CONR[1]*PRIR[1]-PRIR[2];
    double PM;
    if(fabs(coe1) < 1e-15) SM = coe3/coe2;
    else SM = (coe2 - sqrt(coe2*coe2 - 4.*coe1*coe3))/2./coe1;
    PM = (SM*(SL*CONL[2]-CONL[1]) + PRIL[2] - CONL[1]*(SL-PRIL[1])) / (1-SL*SM);
    //std::cout << CONL << "\n" << CONR << std::endl;
    //std::cout << PRIL << "\n" << PRIR << std::endl;
    //std::cout << "coe123" << std::endl;
    //std::cout << coe1 << "\n" << coe2 << "\n" << coe3 << std::endl;
    //std::cout << SM << " " << PM << std::endl;
    bU F;
    F[0] = 0;
    F[1] = PM;
    F[2] = PM*SM;
    //std::cout << "F: " << F << std::endl;
    //double alpha;
    //double tmpl = (1-ul*ul)*csl/(1+fabs(ul)*csl);
    //double tmpr = (1-ur*ur)*csr/(1+fabs(ur)*csr);
    //alpha = std::max(tmpl, tmpr);
    //std::cout << "LF: " << LF(CONL, CONR, PRIL, PRIR, alpha) << std::endl;
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

