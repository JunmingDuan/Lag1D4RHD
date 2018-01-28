#include "Lagranian1D.h"

std::ostream& operator<<(std::ostream& os, const Lagranian1D& H) {
  os.setf(std::ios::scientific);
  os.precision(16);
  for(u_int i=0; i < H.N_x; i++) {
    os << 0.5*(H.mesh[i]+H.mesh[i+1]) << " " << H.Pri[i] << "\n";
  }
  return os;
}

void Lagranian1D::print_con(std::ostream& os) {
  os.setf(std::ios::scientific);
  os.precision(16);
  for(u_int i=0; i < N_x; i++) {
    os << 0.5*(mesh[i]+mesh[i+1]) << " " << Con[i] << "\n";
  }
}

void Lagranian1D::print_pri(std::ostream& os) {
  os.setf(std::ios::scientific);
  os.precision(16);
  for(u_int i=0; i < N_x; i++) {
    os << 0.5*(mesh[i]+mesh[i+1]) << " " << Pri[i] << "\n";
  }
}

void Lagranian1D::print_rupe(std::ostream& os) {
  os.setf(std::ios::scientific);
  os.precision(16);
  for(u_int i=0; i < N_x; i++) {
    os << 0.5*(mesh[i]+mesh[i+1]) << " "
      << Pri[i][0] << " " << Pri[i][1] << " " << Pri[i][2] << " "
      << Pri[i][2]/Pri[i][0]/(Gamma[i]-1) << "\n";
  }
}

double Lagranian1D::fp(const bU& U, double p, const double Gamma) {
  double gamma2 = 1 - pow(U[1]/(p+U[2]),2);
  return U[2]+p - U[0]/sqrt(gamma2) - Gamma/(Gamma-1)*p/gamma2;
}

double Lagranian1D::fpp(const bU& U, double p, const double Gamma) {
  double Ep = U[2] + p;
  double tmp = pow(U[1]/Ep, 2);
  double gamma2 = 1 - tmp;
  return 1 + U[0]*tmp/Ep/pow(gamma2, 1.5) - Gamma/(Gamma-1)*(gamma2 - 2*p*tmp/Ep)/pow(gamma2,2);
}

/*bU Lagranian1D::Con2Pri(const bU& U, const double Gamma) {*/
  //bU prim;
  ////solve a nonlinear equation by Newton method to obtain pressure p
  //u_int ite(0), MAXITE(10);
  //double eps = 1e-15;
  //double a = 0, b = (Gamma-1)*U[2];
  //double p(0), p1(p), y(fp(U,p, Gamma));
  //while(fabs(y) > eps && ite < MAXITE) {
    //p1 = p - y/fpp(U, p, Gamma);
    //y = fp(U, p1, Gamma);
    //ite++;
    //if(fabs(p1-p) < eps) { p = p1; break; }
    //p = p1;
  //}
  //prim[2] = p;
  //prim[1] = U[1]/(U[2]+p);
  //prim[0] = U[0]*sqrt(1-pow(prim[1],2));

  //return prim;
//}

//solve a nonlinear equation by Newton method to obtain pressure p
bU Lagranian1D::Con2Pri(const bU& U, const double Gamma) {
  bU prim;
  u_int ite(0), MAXITE(10);
  double eps = 1e-15;
  double a = 0, b = (Gamma-1)*U[2];
  double p(0.5*b), p1(p), y(fp(U,p, Gamma));
  while(fabs(y) > eps && ite < MAXITE) {
    p1 = p - y/fpp(U, p, Gamma);
    if(fabs(p1-p) < eps) { p = p1; break; }
    if(p1 < a) p1 = a;
    if(p1 > b) p1 = b;
    y = fp(U, p1, Gamma);
    p = p1;
    ite++;
  }
  if(p < 0 || p != p) {
    std::cout << U[0] << " " << U[1] << " " << U[2] << std::endl;
    std::cout << "ite: " << ite << " ; p: " << p << " ; y: " << y << std::endl;
    std::cout << "derivative: " << fpp(U, 0.5*b, Gamma) << std::endl;
    //std::cout << "p1: " << 0.5*b - fp(U,0.5*b,Gamma)/fpp(U, 0.5*b, Gamma) << std::endl;
    std::cout << "f(p1): " << fp(U,p, Gamma) << std::endl;
    std::cout << U[1]/(U[2]+p) << std::endl;
    abort();
  }

  prim[2] = p;
  prim[1] = U[1]/(U[2]+p);
  prim[0] = U[0]*sqrt(1-pow(prim[1],2));
  //std::cout << "ite: " << ite << " ; final p: " << p << std::endl;
  //std::cout << "rho: " << prim[0] << std::endl;
  //abort();

  return prim;
}

bU Lagranian1D::Pri2Con(const bU& U, const double Gamma) {
  bU Con;
  double gamma = 1./sqrt(1-U[1]*U[1]);
  double h = 1 + U[2]/U[0]*Gamma/(Gamma-1);
  Con[0] = gamma*U[0];
  Con[1] = Con[0]*h*gamma*U[1];
  Con[2] = Con[0]*h*gamma - U[2];

  return Con;
}

void Lagranian1D::update_cs(VEC& cs) {
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    cs[i] = cal_cs(Con[i], Pri[i], Gamma[i]);
  }
}

double Lagranian1D::cal_cs
(const bU& Con, const bU& Pri, const double Gamma) {
  double h = 1+Pri[2]/Pri[0]*Gamma/(Gamma-1);
  return sqrt(Gamma*Pri[2]/Pri[0]/h);
}

double Lagranian1D::cal_max_lambda_Lag(int i) {
  double u = Pri[i][1];
  return fabs(1-u*u)*cs[i]/(1+fabs(u)*cs[i]);
}

double Lagranian1D::cal_max_lambda_Eul(int i) {
  double u = Pri[i][1];
  return fabs(1-u*u)*cs[i]/(1+fabs(u)*cs[i])+fabs(u);
}

double Lagranian1D::cal_max_lambda_Lag
(const bU& Con, const bU& Pri, const double Gamma) {
  double h = 1+Pri[2]/Pri[0]*Gamma/(Gamma-1);
  double cs = sqrt(Gamma*Pri[2]/Pri[0]/h);
  double u = Pri[1];
  return fabs(1-u*u)*cs/(1+fabs(u)*cs);
}

double Lagranian1D::cal_max_lambda_Eul
(const bU& Con, const bU& Pri, const double Gamma) {
  double h = 1+Pri[2]/Pri[0]*Gamma/(Gamma-1);
  double cs = sqrt(Gamma*Pri[2]/Pri[0]/h);
  double u = Pri[1];
  return fabs(1-u*u)*cs/(1+fabs(u)*cs)+fabs(cs);
}

double Lagranian1D::cal_max_lambda_Lag() {
  double tmp, a(0);
  for(u_int i = 0; i < N_x; ++i) {
    a = std::max(a, cal_max_lambda_Lag(i));
  }
  return a;
}


bU Lagranian1D::F(const bU& CON, const bU& PRI) {
  bU tmp;
  tmp[0] = 0;
  tmp[1] = PRI[2];
  tmp[2] = PRI[1]*PRI[2];
  return tmp;
}

void Lagranian1D::update_sol(VEC& mesh, Sol& Con, Sol& Pri, Sol& FLUX, const double dt,
    VEC& mesh1, Sol& Con1, Sol& Pri1) {
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con1[i] = (Con[i]*(mesh[i+1]-mesh[i]) - dt*(FLUX[i+1] - FLUX[i])) / (mesh1[i+1]-mesh1[i]);
    Pri1[i] = Con2Pri(Con1[i], Gamma[i]);
  }
}

double Lagranian1D::t_step(const double CFL, double& alpha) {
  double a(1), tmp_lam, hi, tmp_t;
  alpha = 0;
  for(u_int i = 0; i < N_x; ++i) {
    hi = mesh[i+1] - mesh[i];
    tmp_lam = cal_max_lambda_Eul(i);
    alpha = std::max(alpha, tmp_lam);
    tmp_t = hi/alpha;
    a = std::min(a, tmp_t);
  }
  return CFL*a;
}

