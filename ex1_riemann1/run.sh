#!/bin/sh

make -j;
#./main 400 0.4 0.0138858;
./main 400 0.4 0.4;
mv sol.dat ex1_LF_n400_RK2_Lag.dat;

