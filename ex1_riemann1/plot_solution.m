DAT1 = load('ex1_LF_n400_Lag.dat');
x1 = DAT1(:,1);
rho1 = DAT1(:,2);
u1 = DAT1(:,3);
p1 = DAT1(:,4);
e1 = DAT1(:,5);

DAT2 = load('ex1_LF_n400_RK2_Pri_Lag.dat');
x2 = DAT10(:,1);
rho2 = DAT10(:,2);
u2 = DAT10(:,3);
p2 = DAT10(:,4);
e2 = DAT10(:,5);

DAT3 = load('ex1_LF_n400_RK3_Pri_Lag.dat');
x3 = DAT2(:,1);
rho3 = DAT2(:,2);
u3 = DAT2(:,3);
p3 = DAT2(:,4);
e3 = DAT2(:,5);

EX1 = load('../exact_solution/2011GRPex4.2.dat');
x0 = EX1(:,1);
p0 = EX1(:,2);
rho0 = EX1(:,3);
u0 = EX1(:,4);
e0 = EX1(:,5);

plot(x1, rho1/10, '', x2, rho2/10, '', x3, rho3/10, '', x0, rho0/10, '-k');
legend('Order1', 'ENO2', 'WENO3', 'Location', 'NorthEast');
print('ex1_LF_n400_RKC_rho.eps', '-depsc');

plot(x1, u1, '', x2, u2, '', x3, u3, '', x0, u0, '-k');
legend('Order1', 'ENO2', 'WENO3', 'Location', 'NorthEast');
print('ex1_LF_n400_RKC_u.eps', '-depsc');

plot(x1, p1*3/40, '', x2, p2*3/40, '', x3, p3*3/40, '', x0, p0*3/40, '-k');
legend('Order1', 'ENO2', 'WENO3', 'Location', 'NorthEast');
print('ex1_LF_n400_RKC_p.eps', '-depsc');

plot(x1, e1, '', x2, e2, '', x3, e3, '', x0, e0, '-k');
legend('Order1', 'ENO2', 'WENO3', 'Location', 'NorthEast');
print('ex1_LF_n400_RKC_e.eps', '-depsc');


%plot(x3, rho3/10, 'ob', x0, rho0/10, '-k');
%legend('HLLC', 'exact', 'Location', 'NorthEast');
%print('ex1_HLLC_n400_rho.eps', '-depsc');

%plot(x3, u3, 'ob', x0, u0, '-k');
%legend('HLLC', 'exact', 'Location', 'NorthWest');
%print('ex1_HLLC_n400_u.eps', '-depsc');

%plot(x3, p3*3/40, 'ob', x0, p0*3/40, '-k');
%legend('HLLC', 'exact', 'Location', 'NorthEast');
%print('ex1_HLLC_n400_p.eps', '-depsc');

%plot(x3, e3, 'ob', x0, e0, '-k');
%legend('HLLC', 'exact', 'Location', 'NorthEast');
%print('ex1_HLLC_n400_e.eps', '-depsc');

