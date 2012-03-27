function wh_6_21
  w=1; %m/s vertical velocity
  temp=283.15; %K
  mu0=18.27e-6; %Pa s
  T0=291.15; %K
  C=120; %K
  %from http://en.wikipedia.org/wiki/Viscosity
  eta=mu0*(T0 + C)/(temp+C)*(temp/T0)^1.5;  %dynamic viscoisty (Pa s)
  g=9.8; %m/s^2
  rhol=1000.; %kg/m^3
  Vcoeff=2*g*rhol/(9.*eta);
  Gl=100.e-12; %m^2/s
  S=0.005;
  W=1; %m/s
  rcoeff=Gl*S;
  tspan=0:10:1000;
  %
  % redo the 6.21 analytic solution so it starts
  % with a drop of radius 1.e-5 m (10 microns)
  %
  r_analytic=sqrt(2.*Gl*S*tspan + (1.e-5)^2.);
  h_analytic=w*tspan - 2*g*rhol*(Gl*S*tspan.^2. + 1.e-10*tspan)/(9.*eta);
  %start at cloudbase, with a 10 micron drop
  yinit=[0,1.e-5];
  derivs=@(t,y) F(t,y,W,Vcoeff,rcoeff);
  [t,y]=ode45(derivs,tspan,yinit);
  height=y(:,1);
  radius=y(:,2)*1.e6;  %convert radius to microns
  figure(1);
  clf;
  plot(radius,height,'b-','linewidth',3);
  hold on
  plot(r_analytic*1.e6,h_analytic,'g*','markersize',3)
  legend('ode45','analytic');
  xlabel('radius (microns)');
  ylabel('height above cloud base (m)');
  title('WH 6.21: radius vs. height')
  grid on;
end

function yp=F(t,y,W,Vcoeff,rcoeff)
  yp=zeros(2,1); % since output must be a column vector
  %first variable is height in m
  yp(1)=W - Vcoeff*y(2)^2.;
  %second variable is radius in m
  yp(2)= rcoeff/y(2);
end
