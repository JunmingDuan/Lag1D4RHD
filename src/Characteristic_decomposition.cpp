#include "Lagranian1D.h"

void multiply(const bU& x, mvector<mvector<double,3>,3>& M, bU& y) {
  y[0] = M[0][0]*x[0] + M[0][1]*x[1] + M[0][2]*x[2];
  y[1] = M[1][0]*x[0] + M[1][1]*x[1] + M[1][2]*x[2];
  y[2] = M[2][0]*x[0] + M[2][1]*x[1] + M[2][2]*x[2];
}

void Lagranian1D::CHAR_DECOM(const bU& PRIL, const bU& PRIR, const double GAMMAL, const double GAMMAR,
    bU& CHARL, bU& CHARR) {
    double hl, hr, kl, kr, gl, gr, ul, ur;
    double r1, r2, r3;
    hl = 1 + GAMMAL/(GAMMAL-1)*PRIL[2]/PRIL[0];
    hr = 1 + GAMMAR/(GAMMAR-1)*PRIR[2]/PRIR[0];
    kl = sqrt(PRIL[0]*hl);
    kr = sqrt(PRIR[0]*hr);
    ul = PRIL[1];
    ur = PRIR[1];
    gl = 1./sqrt(1-ul*ul);
    gr = 1./sqrt(1-ur*ur);
    r1 = (kl*gl + kr*gr)/(kl+kr);
    r2 = (kl*gl*ul + kr*gr*ur)/(kl+kr);
    r3 = (kl*PRIL[2]/PRIL[0]/hl + kr*PRIR[2]/PRIR[0]/hr) / (kl+kr);
    double cm = 1 - GAMMAL*r3/(GAMMAL-1);
    double cp = 1 + GAMMAL*r3/(GAMMAL-1);
    double v2 = r1*r1 + r2*r2;
    double s2 = 0.5*GAMMAL*r3*(1-v2) - 0.5*(GAMMAL-1)*(1+v2);
    double s = sqrt(s2);
    double e = r1*r1 - r2*r2;
    double y = sqrt((1-GAMMAL*r3)*e + s2);

    mvector<mvector<double,3>,3> M;
    M[0][0] = cm;
    M[0][1] = cm + s2/(GAMMAL-1);
    M[0][2] = cm;
    M[1][0] = r1 - s/y*r2;
    M[1][1] = r1;
    M[1][2] = r1 + s/y*r1;
    M[2][0] = r2 + s/y*r1;
    M[2][1] = r2;
    M[2][2] = r2 - s/y*r2;


}

