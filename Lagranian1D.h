/**
 * @file Lagranian1D.h
 * @brief 1D Lagrangian scheme for relativistic Euler equations
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2018-01-17
 */
#ifndef Lagranian1D_H
#define Lagranian1D_H

#include <iostream>
#include <assert.h>
#include <cstdlib>
#include <omp.h>
#include "para.h"
#include "vvector.h"
#include "mvector.h"

typedef mvector<double,3> bU; //(D, m, E), cell average

class GHOST {
  public:
    bU Con;
    bU Pri;
    double Gamma;
};

class Lagranian1D {
	private:
		typedef bU (*Lfun)(const double t, const double x, double& gamma);
		typedef bU (*NLfun)(bU u, double t, double x);
		typedef vvector<bU> Sol;
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
		Lagranian1D(u_int Nx, double t_start, double t_end, double x_start, double x_end,
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
#pragma omp parallel for num_threads(Nthread)
        for(u_int i = 0; i < N_x+1; ++i) {
          mesh[i] = h*i;
        }
      }

  public:
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

    double fp(const bU& U, double p, const double Gamma);
   double fpp(const bU& U, double p, const double Gamma);
   /**
     * @brief Con2Pri solve (rho,u,p) from (D, m, E)
     *
     * @param U conservative variables (D, m, E)
     * @param gamma
     *
     * @return primitive variables (rho, u, p)
     */
    bU Con2Pri(const bU& U, const double Gamma);
    bU Pri2Con(const bU& U, const double Gamma);
    double cal_max_speed(int i);

    double t_step(const double CFL, double& alpha) {
      double a(1), tmp_lam, hi, tmp_t;
      alpha = 0;
      for(u_int i = 0; i < N_x; ++i) {
        hi = mesh[i+1] - mesh[i];
        tmp_lam = cal_max_speed(i);
        if(tmp_lam > alpha) alpha = tmp_lam;
        tmp_t = hi/alpha;
        if(tmp_t < a) a = tmp_t;
      }
      return CFL*a;
    }

    bU F(const bU& CON, const bU& PRI);

    bU LF(const bU& CONL, const bU& CONR, const bU& PRIL, const bU& PRIR, const double alpha);
    void cal_flux_LF(double alpha);
    void forward_LF(double dt, double alpha);

    bU HLLC(const bU& CONL, const bU& CONR, const bU& PRIL, const bU& PRIR, const double Gammal, const double Gammar);
    void cal_flux_HLLC(double alpha);
    void forward_HLLC(double dt, double alpha);

    void cal_us_roeav();
    void move_mesh(double dt);

  public:
		void Solve() {
      double t_now(t_start), dt(0), alpha(0);
			Initial();

      int n = 0;
      while(t_now < t_end - 1e-10) {
      //while( n < 1 ) {
        dt = t_step(CFL, alpha);
        if(dt + t_now > t_end) dt = t_end - t_now;

        LF(Con[0], Con[1], Pri[0], Pri[1], 1);
        //forward_LF(dt, alpha);
        //forward_HLLC(dt, alpha);

        t_now += dt;
        std::cout << "t: " << t_now << " , dt: " << dt << std::endl;
        n ++ ;
      }
		}

    void print_con(std::ostream& os);
    void print_pri(std::ostream& os);
    void print_rupe(std::ostream& os);
		friend std::ostream& operator<<(std::ostream&, const Lagranian1D&);
};

#endif

