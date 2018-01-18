/**
 * @file Lagranian1D.h
 * @brief 1D Lagrangian scheme for relativistic Euler equations
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2018-01-17
 */
#ifndef SCL1D_H
#define SCL1D_H

#include <iostream>
#include <assert.h>
#include <cstdlib>
#include "vvector.h"

typedef mvector<double,3> bU; //(D, m, E), cell average

double fp(const bU& U, double p, const double Gamma) {
  double gamma2 = 1 - pow(U[1]/(p+U[2]),2);
  return U[2]+p - U[0]/sqrt(gamma2) - Gamma/(Gamma-1)*p/gamma2;
}
double fpp(const bU& U, double p, const double Gamma) {
  double Ep = U[2] + p;
  double tmp = pow(U[1]/Ep, 2);
  double gamma2 = 1 - tmp;
  return 1 + U[0]*tmp/Ep/pow(gamma2, 1.5) - Gamma/(Gamma-1)*(gamma2 - 2*p*tmp/Ep)/pow(gamma2,2);
}
/**
 * @brief Con2Pri solve (rho,u,p) from (D, m, E)
 *
 * @param U conservative variables (D, m, E)
 * @param gamma
 *
 * @return primitive variables (rho, u, p)
 */
bU Con2Pri(const bU& U, const double Gamma) {
  bU prim;
  //solve a nonlinear equation by Newton method to obtain pressure p
  u_int ite(0), MAXITE(20);
  double eps = 1e-13;
  double p(0), p1(p), y(fp(U,p, Gamma));
  while(fabs(y) > eps && ite < MAXITE) {
    p1 = p - y/fpp(U, p, Gamma);
    y = fp(U, p1, Gamma);
    ite++;
    if(fabs(p1-p) < eps) { p = p1; break; }
    p = p1;
  }
  //std::cout << ite << " " << y << std::endl;
  prim[2] = p;
  prim[1] = U[1]/(U[2]+p);
  prim[0] = U[0]*sqrt(1-pow(prim[1],2));

  return prim;
}

bU Pri2Con(const bU& U, const double Gamma) {
  bU Con;
  double gamma = 1./sqrt(1-U[1]*U[1]);
  double h = 1 + U[2]/U[0]*Gamma/(Gamma-1);
  Con[0] = gamma*U[0];
  Con[1] = Con[0]*h*gamma*U[1];
  Con[2] = Con[0]*h*gamma - U[2];

  return Con;
}

class GHOST {
  public:
    bU Con;
    bU Pri;
    double Gamma;
};

template<class T>
class SCL1D {
	private:
		typedef T (*Lfun)(const double t, const double x, double& gamma);
		typedef T (*NLfun)(T u, double t, double x);
		typedef vvector<T> Sol;
		u_int N_x;
		double t_start;
		double t_end;
		double x_start;
		double x_end;
		Lfun initial;
		double CFL;
		Sol Con;//数值解,conservative variables
		Sol Pri;//数值解,primitive variables
    vvector<double> Di;
    vvector<double> mesh;
    vvector<double> Gamma;
    vvector<double> us;
    Sol FLUX;
    GHOST ghostl, ghostr;

	public:
		SCL1D(u_int Nx, double t_start, double t_end, double x_start, double x_end,
        Lfun initial, double CFL = 0.5)
      : N_x(Nx), t_start(t_start), t_end(t_end), x_start(x_start), x_end(x_end),
      initial(initial), CFL(CFL) {
        Con.resize(N_x);
        Pri.resize(N_x);
        Di.resize(N_x);
        mesh.assign(N_x+1, 0);
        Gamma.assign(N_x, 0);
        us.assign(N_x+1, 0);
        FLUX.resize(N_x+1);
        double h = (x_end - x_start) / N_x;
        for(u_int i = 0; i < N_x+1; ++i) {
          mesh[i] = h*i;
        }
      }

  private:
    void Initial() {
      for(u_int j = 0; j != N_x; ++j) {
        Pri[j] = initial(t_start, 0.5*(mesh[j]+mesh[j+1]), Gamma[j]);
        Con[j] = Pri2Con(Pri[j], Gamma[j]);
      }
    }

    void InfiniteBD() {
      ghostl.Con = Con[0];
      ghostl.Pri = Pri[0];
      ghostr.Con = Con[N_x-1];
      ghostr.Pri = Pri[N_x-1];
      ghostl.Gamma = Gamma[0];
      ghostr.Gamma = Gamma[N_x-1];
    }

    double cal_speed(int i) {
      double h = 1+Pri[i][2]/Pri[i][0]*Gamma[i]/(Gamma[i]-1);
      double cs = sqrt(Gamma[i]*Pri[i][2]/Pri[i][0]/h);
      double u = Pri[i][1];
      double lam1 = (u-cs)/(1-u*cs)-u;
      double lam3 = (u+cs)/(1+u*cs)-u;
      return std::max(fabs(lam1), fabs(lam3));
    }

    double t_step(const double CFL, double& alpha) {
      double a(1e-2), tmp_lam, hi, tmp_t;
      alpha = 0;
      for(u_int i = 0; i < N_x; ++i) {
        hi = mesh[i+1] - mesh[i];
        tmp_lam = cal_speed(i);
        if(tmp_lam > alpha) alpha = tmp_lam;
        tmp_t = hi/alpha;
        if(tmp_t < a) a = tmp_t;
      }
      //std::cout << "alpha: " << alpha << std::endl;
      //std::cout << "a: " << a << std::endl;
      return CFL*a;
    }

    //void rp_solver(const bU& UL, const bU& UR, double& us, double& ps, const double gammal, const double gammar) {
      //double pl = EOS(UL, gammal);
      //double pr = EOS(UR, gammar);
      //double ul = UL[1];
      //double ur = UR[1];
      //double al = sqrt(gammal*pl/UL[0]);
      //double ar = sqrt(gammar*pr/UR[0]);
      //if(al != al) std::cout << "al wrong: " << pl << " " << UL[0] << std::endl;
      //if(ar != ar) std::cout << "ar wrong: " << pr << " " << UR[0] << std::endl;
      //double zl = UL[0] * al;
      //double zr = UR[0] * ar;
      //double sum = zl + zr;
      //us = (zl*ul + zr*ur)/sum - (pr-pl)/sum;
      //ps = (zl*pr + zr*pl)/sum - zl*zr*(ur-ul)/sum;
    //}

    bU F(const bU& CON, const bU& PRI) {
      bU tmp;
      tmp[0] = 0;
      tmp[1] = PRI[2];
      tmp[2] = PRI[1]*PRI[2];
      return tmp;
    }

    bU LF(const bU& CONL, const bU& CONR, const bU& PRIL, const bU& PRIR, const double alpha) {
      return 0.5*(F(CONL, PRIL) + F(CONR, PRIR)) - 0.5*(CONR-CONL)*alpha;
    }

    void cal_flux(double alpha) {
      for(u_int i = 0; i < N_x+1; ++i) {
        if(i == 0) FLUX[i] = LF(ghostl.Con, Con[i], ghostl.Pri, Pri[i], alpha);
        else if(i == N_x) FLUX[i] = LF(Con[i-1], ghostr.Con, Pri[i-1], ghostr.Pri, alpha);
        else FLUX[i] = LF(Con[i-1], Con[i], Pri[i-1], Pri[i], alpha);
      }
    }

    void cal_us() {
      double srhol, srhor, ul, ur;
      for(u_int i = 0; i < N_x+1; ++i) {
        if(i == 0) {
          srhol = sqrt(ghostl.Pri[0]);
          ul = sqrt(ghostl.Pri[1]);
        }
        else {
          srhol = sqrt(Pri[i-1][0]);
          ul = Pri[i-1][1];
        }
        if(i == N_x) {
          srhor = sqrt(ghostr.Pri[0]);
          ur = sqrt(ghostr.Pri[1]);
        }
        else {
          srhor = sqrt(Pri[i][0]);
          ur = Pri[i][1];
        }
        us[i] = (srhol*ul + srhor*ur)/(srhol+srhor);
      }
      us[0] = 0;
      us[N_x] = 0;
    }

    void move_mesh(double dt) {
      for(u_int i = 0; i < N_x+1; ++i) {
        mesh[i] += dt * us[i];
      }
    }

    void forward(double dt, double alpha) {
      InfiniteBD();
      cal_flux(alpha);
      cal_us();
      for(u_int i = 0; i < N_x; ++i) {
        Con[i] *= (mesh[i+1]-mesh[i]);
      }
      move_mesh(dt);
      for(u_int i = 0; i < N_x; ++i) {
        Con[i] = (Con[i] - dt*(FLUX[i+1] - FLUX[i])) / (mesh[i+1]-mesh[i]);
      }
      for(u_int i = 0; i < N_x; ++i) {
        Pri[i] = Con2Pri(Con[i], Gamma[i]);
      }
    }

  public:
		void Solve() {
      double t_now(t_start), dt(0), alpha(0);
			Initial();

      while(t_now < t_end - 1e-10) {
        dt = t_step(CFL, alpha);
        if(dt + t_now > t_end) dt = t_end - t_now;

        forward(dt, alpha);

        t_now += dt;
        std::cout << "t: " << t_now << " , dt: " << dt << std::endl;
      }
		}

		template<class Type> friend std::ostream& operator<<(std::ostream&,SCL1D<Type>&);

    void print_con(std::ostream& os) {
      os.setf(std::ios::scientific);
      os.precision(16);
      for(u_int i=0; i < N_x; i++) {
        os << 0.5*(mesh[i]+mesh[i+1]) << " " << Con[i] << "\n";
      }
    }

    void print_pri(std::ostream& os) {
      os.setf(std::ios::scientific);
      os.precision(16);
      for(u_int i=0; i < N_x; i++) {
        os << 0.5*(mesh[i]+mesh[i+1]) << " " << Pri[i] << "\n";
      }
    }

    void print_rupe(std::ostream& os) {
      os.setf(std::ios::scientific);
      os.precision(16);
      for(u_int i=0; i < N_x; i++) {
        os << 0.5*(mesh[i]+mesh[i+1]) << " "
          << Pri[i][0] << " " << Pri[i][1] << " " << Pri[i][2] << " "
          << Pri[i][2]/Pri[i][0]/(Gamma[i]-1) << "\n";
      }
    }

};

template<class T>
std::ostream& operator<<(std::ostream& os, SCL1D<T>& H){
	os.setf(std::ios::scientific);
  os.precision(16);
	for(u_int i=0; i < H.N_x; i++) {
		os << 0.5*(H.mesh[i]+H.mesh[i+1]) << " " << H.Pri[i] << "\n";
	}
	return os;
}

#endif

//std::cout << "conservative" << std::endl;
//print_con(std::cout);
//std::cout << "primitive" << std::endl;
//print_pri(std::cout);
//for(u_int i=0; i < N_x; i++) {
//std::cout << 0.5*(mesh[i]+mesh[i+1]) << " " << Con2Pri(Con[i], Gamma[i]) << "\n";
//}


