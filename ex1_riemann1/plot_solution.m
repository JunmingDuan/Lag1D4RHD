%DAT1 = load('sol.dat');
DAT1 = load('ex1_LF_n400_Lag.dat');
x1 = DAT1(:,1);
rho1 = DAT1(:,2);
u1 = DAT1(:,3);
p1 = DAT1(:,4);
e1 = DAT1(:,5);

DAT10 = load('ex1_LF_n400_Eul.dat');
x10 = DAT10(:,1);
rho10 = DAT10(:,2);
u10 = DAT10(:,3);
p10 = DAT10(:,4);
e10 = DAT10(:,5);

DAT2 = load('ex1_LLF_n400.dat');
x2 = DAT2(:,1);
rho2 = DAT2(:,2);
u2 = DAT2(:,3);
p2 = DAT2(:,4);
e2 = DAT2(:,5);

DAT20 = load('ex1_LLF_n400_Eul.dat');
x20 = DAT20(:,1);
rho20 = DAT20(:,2);
u20 = DAT20(:,3);
p20 = DAT20(:,4);
e20 = DAT20(:,5);

DAT3 = load('ex1_HLLC_n400.dat');
x3 = DAT3(:,1);
rho3 = DAT3(:,2);
u3 = DAT3(:,3);
p3 = DAT3(:,4);
e3 = DAT3(:,5);

EX1 = load('../exact_solution/2011GRPex4.2.dat');
x0 = EX1(:,1);
p0 = EX1(:,2);
rho0 = EX1(:,3);
u0 = EX1(:,4);
e0 = EX1(:,5);

%{%ex1%}
%figure(1)
%plot(x1, rho1/10, 'or', x2, rho2/10, '-b', x3, rho3/10, '*b', x0, rho0/10, '-k');
%%legend('Recon', 'unRecon', 'exact');
%legend('LF', 'LLF', 'HLLC', 'exact');
%figure(2)
%plot(x1, u1, 'or', x2, u2, '-b', x3, u3, '*b', x0, u0, '-k');
%%legend('Recon', 'unRecon', 'exact');
%legend('LF', 'LLF', 'HLLC', 'exact');
%figure(3)
%plot(x1, p1*3/40, 'or', x2, p2*3/40, '-b', x3, p3*3/40, '*b', x0, p0*3/40, '-k');
%%legend('Recon', 'unRecon', 'exact');
%legend('LF', 'LLF', 'HLLC', 'exact');
%figure(4)
%plot(x1, e1, 'or', x2, e2, '-b', x3, e3, '*b', x0, e0, '-k');
%%legend('Recon', 'unRecon', 'exact');
%legend('LF', 'LLF', 'HLLC', 'exact');

%ex1
%plot(x1, rho1/10, 'or', x10, rho10/10, '*b', x0, rho0/10, '-k');
%legend('Lag', 'EUl', 'Location', 'NorthEast');
%print('ex1_LF_n400_rho.eps', '-depsc');

%plot(x1, u1, 'or', x10, u10, '*b', x0, u0, '-k');
%legend('Lag', 'EUl', 'Location', 'NorthWest');
%print('ex1_LF_n400_u.eps', '-depsc');

%plot(x1, p1*3/40, 'or', x10, p10*3/40, '*b', x0, p0*3/40, '-k');
%legend('Lag', 'EUl', 'Location', 'NorthEast');
%print('ex1_LF_n400_p.eps', '-depsc');

%plot(x1, e1, 'or', x10, e10, '*b', x0, e0, '-k');
%legend('Lag', 'EUl', 'Location', 'NorthEast');
%print('ex1_LF_n400_e.eps', '-depsc');

plot(x2, rho2/10, 'or', x20, rho20/10, '*b', x0, rho0/10, '-k');
legend('Lag', 'Eul', 'Location', 'NorthEast');
print('ex1_LLF_n400_rho.eps', '-depsc');

plot(x2, u2, 'or', x20, u20, '*b', x0, u0, '-k');
legend('Lag', 'Eul', 'Location', 'NorthWest');
print('ex1_LLF_n400_u.eps', '-depsc');

plot(x2, p2*3/40, 'or', x20, p20*3/40, '*b', x0, p0*3/40, '-k');
legend('Lag', 'Eul', 'Location', 'NorthEast');
print('ex1_LLF_n400_p.eps', '-depsc');

plot(x2, e2, 'or', x20, e20, '*b', x0, e0, '-k');
legend('Lag', 'Eul', 'Location', 'NorthEast');
print('ex1_LLF_n400_e.eps', '-depsc');


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

