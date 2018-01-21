#include "Lagranian1D.h"

void Lagranian1D::Euler_forward_LF(double dt, double alpha) {
  InfiniteBD(Con, Pri);
  cal_flux_LF(Con, Pri, FLUX, alpha);
  cal_us_roeav(Con, Pri, us);
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] *= (mesh[i+1]-mesh[i]);
  }
  move_mesh(mesh, us, dt, mesh);
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] = (Con[i] - dt*(FLUX[i+1] - FLUX[i])) / (mesh[i+1]-mesh[i]);
  }
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Pri[i] = Con2Pri(Con[i], Gamma[i]);
  }
}

void Lagranian1D::Euler_forward_HLLC(const double dt) {
  InfiniteBD(Con, Pri);
  cal_flux_HLLC(Con, Pri, FLUX);
  cal_us_roeav(Con, Pri, us);
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] *= (mesh[i+1]-mesh[i]);
  }
  move_mesh(mesh, us, dt, mesh);
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

void Lagranian1D::SSP_RK_LF(Sol& Con, Sol& Pri, vvector<double>& mesh, const double dt, double alpha) {
  vvector<double> mesh_n(mesh), mesh1(N_x+1);
  Sol Con_n(Con), Pri_n(Pri);
  //stage 1
  InfiniteBD(Con_n, Pri_n);
  cal_flux_LF(Con_n, Pri_n, FLUX, alpha);
  cal_us_roeav(Con_n, Pri_n, us);
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x+1; ++i) {
    mesh[i] = mesh[i] + dt * us[i];
  }
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] = (Con_n[i]*(mesh_n[i+1]-mesh_n[i]) - dt*(FLUX[i+1] - FLUX[i])) / (mesh[i+1]-mesh[i]);
    Pri[i] = Con2Pri(Con[i], Gamma[i]);
  }
  //stage 2
  InfiniteBD(Con, Pri);
  cal_flux_LF(Con, Pri, FLUX, alpha);
  cal_us_roeav(Con, Pri, us);
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x+1; ++i) {
    mesh1[i] = 0.75*mesh_n[i] + 0.25*(mesh[i] + dt * us[i]);
  }
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] = ( 0.75*Con_n[i]*(mesh_n[i+1]-mesh_n[i])
      + 0.25*(Con[i]*(mesh[i+1]-mesh[i]) - dt*(FLUX[i+1] - FLUX[i])) )/ (mesh1[i+1]-mesh1[i]);
    Pri[i] = Con2Pri(Con[i], Gamma[i]);
  }
  //stage 3
  InfiniteBD(Con, Pri);
  cal_flux_LF(Con, Pri, FLUX, alpha);
  cal_us_roeav(Con, Pri, us);
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x+1; ++i) {
    mesh[i] = 1./3*mesh_n[i] + 2./3*(mesh1[i] + dt * us[i]);
  }
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] = ( 1./3*Con_n[i]*(mesh_n[i+1]-mesh_n[i])
      + 2./3*(Con[i]*(mesh1[i+1]-mesh1[i]) - dt*(FLUX[i+1] - FLUX[i])) )/ (mesh[i+1]-mesh[i]);
    Pri[i] = Con2Pri(Con[i], Gamma[i]);
  }

}

void Lagranian1D::SSP_RK_HLLC(Sol& Con, Sol& Pri, vvector<double>& mesh, const double dt) {
  vvector<double> mesh_n(mesh);
  Sol Con_n = Con, Pri_n = Pri;
  //stage 1
  InfiniteBD(Con_n, Pri_n);
  cal_flux_HLLC(Con_n, Pri_n, FLUX);
  cal_us_roeav(Con_n, Pri_n, us);
  move_mesh(mesh_n, us, dt, mesh);
  update_sol(mesh_n, Con_n, Pri_n, FLUX, dt, mesh, Con, Pri);
  //stage 2
}

