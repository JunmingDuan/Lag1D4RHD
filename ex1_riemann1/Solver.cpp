#include "Lagranian1D.h"

void Lagranian1D::Solver() {
  double t_now(t_start), dt(0), alpha(0);
  Initial();
  update_cs(cs);

  int n = 0;
  while(t_now < t_end) {
    //while( n < 1 ) {
    dt = t_step(CFL, alpha);
    dt = std::min(dt, t_end-t_now);
    if(dt != dt) abort();

    Euler_forward_LF(dt, alpha, mesh);
    //Euler_forward_LLF(dt, mesh);
    //Euler_forward_HLLC(dt, mesh);
    //RK2_LF(Con, Pri, mesh, dt, alpha);
    //SSP_RK_LF(Con, Pri, mesh, dt, alpha);
    //SSP_RK_HLLC(Con, Pri, mesh, dt);

    t_now += dt;
    std::cout << "t: " << t_now << " , dt: " << dt << std::endl;
    n ++ ;
  }
}

