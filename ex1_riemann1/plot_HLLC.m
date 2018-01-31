%DAT1 = load('ex1_HLLC_n400_RK2_Cha_Lag.dat');
DAT1 = load('ex1_HLLC_n400_RK3_Con_Lag.dat');
x1 = DAT1(:,1);
rho1 = DAT1(:,2);
u1 = DAT1(:,3);
p1 = DAT1(:,4);
e1 = DAT1(:,5);

EX1 = load('../exact_solution/2011GRPex4.2.dat');
x0 = EX1(:,1);
p0 = EX1(:,2);
rho0 = EX1(:,3);
u0 = EX1(:,4);
e0 = EX1(:,5);

figure(1)
plot(x1, rho1/10, 'or', x0, rho0/10, '-k');
legend('Lag', 'Eul', 'Location', 'NorthEast');
figure(2)
plot(x1, u1, 'or',  x0, u0, '-k');
legend('Lag', 'Eul', 'Location', 'NorthWest');
figure(3)
plot(x1, p1*3/40, 'or',  x0, p0*3/40, '-k');
legend('Lag', 'Eul', 'Location', 'NorthEast');
figure(4)
plot(x1, e1, 'or', x0, e0, '-k');
legend('Lag', 'Eul', 'Location', 'NorthEast');

%figure(1)
%plot(x1, rho1/10, 'or', x10, rho10/10, '*b', x0, rho0/10, '-k');
%legend('Lag', 'Eul', 'Location', 'NorthEast');
%figure(2)
%plot(x1, u1, 'or', x10, u10, '*b', x0, u0, '-k');
%legend('Lag', 'Eul', 'Location', 'NorthWest');
%figure(3)
%plot(x1, p1*3/40, 'or', x10, p10*3/40, '*b', x0, p0*3/40, '-k');
%legend('Lag', 'Eul', 'Location', 'NorthEast');
%figure(4)
%plot(x1, e1, 'or', x10, e10, '*b', x0, e0, '-k');
%legend('Lag', 'Eul', 'Location', 'NorthEast');

