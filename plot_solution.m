DAT = load('sol.dat');
x = DAT(:,1);
rho = DAT(:,2);
u = DAT(:,3);
p = DAT(:,4);
e = DAT(:,5);

DAT1 = load('LF.dat');
x1 = DAT1(:,1);
rho1 = DAT1(:,2);
u1 = DAT1(:,3);
p1 = DAT1(:,4);
e1 = DAT1(:,5);

%ex1
figure(1)
plot(x, rho/10, '-o', x1, rho1/10, '-');
figure(2)
plot(x, u, '-o', x1, u1, '-');
figure(3)
plot(x, p*3/40, '-o', x1, p1*3/40, '-');
figure(4)
plot(x, e, '-o', x1, e1, '-');


%ex3
%figure(1)
%plot(x, rho, '-o');
%axis([0.49,0.54,-10,130]);
%figure(2)
%plot(x, u, '-o');
%axis([0.49,0.54,-1,1.2]);
%figure(3)
%plot(x, p, '-o');
%axis([0.49,0.54,-25,400]);
%figure(4)
%plot(x, e, '-o');
%axis([0.49,0.54,-50,800]);


