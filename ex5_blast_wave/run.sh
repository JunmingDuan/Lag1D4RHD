#!/bin/sh

make -j;
./main 1000 0.4 0.43;
#./main 400 0.4 0.411;
#mv sol.dat ex5_LLF_n400_RK1_Lag.dat;
#mv sol.dat ex5_LF_n400_RK3_Cha_Lag.dat;
#mv sol.dat ex5_HLLC_n400_RK3_Cha_Lag.dat;

