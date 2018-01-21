#include "Lagranian1D.h"

void Lagranian1D::Solver() {
  double t_now(t_start), dt(0), alpha(0);
  Initial();

  int n = 0;
  while(t_now < t_end) {
    //while( n < 1 ) {
    update_cs(cs);
    dt = t_step(CFL, alpha);
    dt = std::min(dt, t_end-t_now);

    //Euler_forward_LF(dt, alpha);
    //Euler_forward_LLF(dt);
    Euler_forward_HLLC(dt);
    //SSP_RK_LF(Con, Pri, mesh, dt, alpha);
    //SSP_RK_HLLC(Con, Pri, mesh, dt);

    t_now += dt;
    std::cout << "t: " << t_now << " , dt: " << dt << std::endl;
    n ++ ;
  }
}

