#!/bin/sh

make -j;
./main 400 0.4 0.4;
#mv sol.dat ex1_LLF_n400_RK1_Lag.dat;
#mv sol.dat ex1_LF_n400_RK3_Cha_Lag.dat;
#mv sol.dat ex1_HLLC_n400_RK3_Cha_Lag.dat;

