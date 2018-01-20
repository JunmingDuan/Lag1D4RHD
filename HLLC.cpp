#include "Lagranian1D.h"

bU Lagranian1D::HLLC(const bU& CONL, const bU& CONR, const bU& PRIL, const bU& PRIR, const double Gammal, const double Gammar, double& SM) {
  double SL, SR, ui, ci, vl, vr;
  double srhol, srhor;
  double hl, csl, ul, hr, csr, ur;
  double lam1, lam2, lam3, lam4;

  //characteristic speed of left and right side
  hl = 1+PRIL[2]/PRIL[0]*Gammal/(Gammal-1);
  csl = sqrt(Gammal*PRIL[2]/PRIL[0]/hl);
  //std::cout << "hl,p,rho: " << hl << " " << PRIL[2] << " " << PRIL[0] << std::endl;
  ul = PRIL[1];
  lam1 = -(1-ul*ul)*csl/(1-ul*csl);//negative speed
  lam2 = (1-ul*ul)*csl/(1+ul*csl);//positive speed
  hr = 1+PRIR[2]/PRIR[0]*Gammar/(Gammar-1);
  csr = sqrt(Gammar*PRIR[2]/PRIR[0]/hr);
  ur = PRIR[1];
  lam3 = -(1-ur*ur)*csr/(1-ur*csr);//negative speed
  lam4 = (1-ur*ur)*csr/(1+ur*csr);//positive speed
  SL = std::min((lam1), (lam3));//left characteristic speed
  SR = std::max((lam2), (lam4));//right characteristic speed
  //roe_average of ul, ur
  //srhol = sqrt(PRIL[0]);
  //srhor = sqrt(PRIR[0]);
  //ui = (srhol*ul + srhor*ur)/(srhol+srhor);
  //left and right speed limit
  //vl = ul - ui;
  //vr = ur - ui;

  //S1 = std::min(vl-cl, -ci);
  //S2 = std::max(vr+cr, ci);

  //std::cout << "left and right speed: " << SL << " " << SR << std::endl;
  if(SL > 0) { return F(CONL, PRIL); }
  else if(SR < 0) {  return F(CONR, PRIR); }
  else {
    //cal S^* and p^*, S^* is the smaller root of a quadric equation:
    //F_E^{hll} x^2 - (E^{hll}+F^{hll}_m) x + m^{hll} = 0
    double coe1 = SR*CONL[1] - SL*CONR[1] + SL*SR*(CONR[2] - CONL[2]);
    double coe2 = ( SR*(CONL[1]*PRIL[1]+PRIL[2]) - SL*(CONR[1]*PRIR[1]+PRIR[2]) + SL*SR*(CONR[1] - CONL[1])
        + SR*CONR[2] - SL*CONL[2] + CONL[1] - CONR[1] );
    double coe3 = SR*CONR[1] - SL*CONL[1] + CONL[1]*PRIL[1]+PRIL[2] - CONR[1]*PRIR[1]-PRIR[2];
    //std::cout << "SL,SR: " << SL << " " << SR << std::endl;
    //std::cout << "COE: " << coe1 << " " << coe2 << " " << coe3 << std::endl;
    double PM;
    if(fabs(coe1) < 1e-10) SM = coe3/coe2;
    else SM = (coe2 - sqrt(coe2*coe2 - 4.*coe1*coe3))/2./coe1;
    PM = (SM*(SL*CONL[2]-CONL[1]) + PRIL[2] - CONL[1]*(SL-PRIL[1])) / (1+SL*SM);
    //std::cout << "star wave" << std::endl;
    bU F;
    F[0] = 0;
    F[1] = PM;
    F[2] = PM*SM;
    //std::cout << "PM: " << PM << " ; SM: " << SM << std::endl;
    return F;
  }
}

void Lagranian1D::cal_flux_HLLC() {
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x+1; ++i) {
    if(i == 0) FLUX[i] = HLLC(ghostl.Con, Con[i], ghostl.Pri, Pri[i], ghostl.Gamma, Gamma[i], us[i]);
    else if(i == N_x) FLUX[i] = HLLC(Con[i-1], ghostr.Con, Pri[i-1], ghostr.Pri, Gamma[i-1], ghostr.Gamma, us[i]);
    else {
      FLUX[i] = HLLC(Con[i-1], Con[i], Pri[i-1], Pri[i], Gamma[i-1], Gamma[i], us[i]);
    }
  }
}

void Lagranian1D::forward_HLLC(double dt, double alpha) {
  InfiniteBD();
  cal_us_roeav();
  cal_flux_HLLC();
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] *= (mesh[i+1]-mesh[i]);
  }
  move_mesh(dt);
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] = (Con[i] - dt*(FLUX[i+1] - FLUX[i])) / (mesh[i+1]-mesh[i]);
    if(Con[i][0] != Con[i][0]) {
      std::cout << "NAN, D: " << Con[i][0] << std::endl;
      std::cout << "i: " << i << std::endl;
      abort();
    }
    if(Con[i][1] != Con[i][1]) {
      std::cout << "NAN, m: " << Con[i][1] << std::endl;
      std::cout << "i: " << i << std::endl;
      abort();
    }
    if(Con[i][2] != Con[i][2]) {
      std::cout << "NAN, E: " << Con[i][2] << std::endl;
      std::cout << "i: " << i << std::endl;
      abort();
    }
  }
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Pri[i] = Con2Pri(Con[i], Gamma[i]);
    if(Pri[i][2] != Pri[i][2]) {
      std::cout << "NAN, p: " << Pri[i][2] << std::endl;
      std::cout << "i: " << i << std::endl;
      abort();
    }
    if(Pri[i][1] != Pri[i][1]) {
      std::cout << "NAN, u: " << Pri[i][1] << std::endl;
      std::cout << "i: " << i << std::endl;
      abort();
    }
    if(Pri[i][0] != Pri[i][0]) {
      std::cout << "NAN, r: " << Pri[i][0] << std::endl;
      std::cout << "i: " << i << std::endl;
      std::cout << "u: " << Pri[i][1] << std::endl;
      abort();
    }
  }
}


