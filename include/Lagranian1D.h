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
    vvector<double> cs;
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
        cs.assign(N_x, 0);
        FLUX.resize(N_x+1);
        double h0 = (x_end - x_start) / N_x;
#pragma omp parallel for num_threads(Nthread)
        for(u_int i = 0; i < N_x+1; ++i) {
          mesh[i] = h0*i;
        }
      }

  public:
    void Initial() {
      for(u_int j = 0; j != N_x; ++j) {
        Pri[j] = initial(t_start, 0.5*(mesh[j]+mesh[j+1]), Gamma[j]);
        Con[j] = Pri2Con(Pri[j], Gamma[j]);
      }
    }

    void InfiniteBD(Sol& Con, Sol& Pri) {
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
    double cal_cs(const bU& Con, const bU& Pri, const double Gamma);
    void update_cs(vvector<double>&);
    double cal_max_lambda_Lag(int i);
    double cal_max_lambda_Eul(int i);
    double cal_max_lambda_Lag(const bU& Con, const bU& Pri, const double Gamma);
    double cal_max_lambda_Eul(const bU& Con, const bU& Pri, const double Gamma);

    double t_step(const double CFL, double& alpha);

    bU F(const bU& CON, const bU& PRI);

    bU LF(const bU& CONL, const bU& CONR, const bU& PRIL, const bU& PRIR, const double alpha);
    void cal_flux_LF(Sol& Con, Sol& Pri, Sol& FLUX, double alpha);

    bU LLF(const bU& CONL, const bU& CONR, const bU& PRIL, const bU& PRIR, const double alpha);
    void cal_flux_LLF(Sol& Con, Sol& Pri, Sol& FLUX);

    bU HLLC(const bU& CONL, const bU& CONR, const bU& PRIL, const bU& PRIR, const double Gammal, const double Gammar, double&);
    void cal_flux_HLLC(Sol& Con, Sol& Pri, Sol& FLUX, vvector<double>& us);

    void cal_us_roeav(Sol& Con, Sol& Pri, vvector<double>& us);
    void move_mesh(vvector<double>&, vvector<double>&, double dt, vvector<double>&);

    void Euler_forward_LF(double dt, double alpha);
    void Euler_forward_LLF(double dt);
    void Euler_forward_HLLC(const double dt);

    void SSP_RK_LF(Sol& Con, Sol& Pri, vvector<double>& mesh, const double dt, double alpha);
    void SSP_RK_HLLC(Sol& Con, Sol& Pri, vvector<double>& mesh, const double dt);

  public:
    void Solver();

    void update_sol(vvector<double>& mesh, Sol& Con, Sol& Pri, Sol& FLUX, const double dt,
        vvector<double>& mesh1, Sol& Con1, Sol& Pri1);

    void print_con(std::ostream& os);
    void print_pri(std::ostream& os);
    void print_rupe(std::ostream& os);
    friend std::ostream& operator<<(std::ostream&, const Lagranian1D&);
};

#endif

