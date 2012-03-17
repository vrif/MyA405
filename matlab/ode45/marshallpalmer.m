D=0:.1:5;
Dmm=D;
Dcm=D*0.1;
N0=0.08*1.e6*1.e-1;%m^{-3} mm^{-1}
R=1.;
theLambda=41.*R^(-0.21);
curve1=N0*exp(-theLambda*Dcm);
R=5.;
theLambda=41.*R^(-0.21);
curve2=N0*exp(-theLambda*Dcm);
R=25.;
theLambda=41.*R^(-0.21);
curve3=N0*exp(-theLambda*Dcm);
figure(1);
semilogy(D,curve1);
hold on;
semilogy(D,curve2,'r-');
semilogy(D,curve3,'g-');
xlabel('Drop diameter (mm)')
ylabel('n(D) m^{-3} mm^{-1}')
title('Marshall Palmer distribution for three rain rates');
legend(' 1 mm/hr','5 mm/hr','25 mm/hr')
print -depsc marsallpalmer.eps
hold off;
figure(2);
plot(D,curve1);
hold on;
plot(D,curve2,'r-');
plot(D,curve3,'g-');
xlabel('Drop diameter (mm)')
ylabel('n(D) m^{-3} mm^{-1}')
title('Marshall Palmer distribution for three rain rates');
legend(' 1 mm/hr','5 mm/hr','25 mm/hr')
print -depsc marsallpalmer.eps
hold off;
