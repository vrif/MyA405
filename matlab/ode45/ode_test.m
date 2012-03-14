function ode_test
  tspan=0:0.01:10;
  [t,y]=ode45(@F,tspan,[2,8]);
  figure(1);
  clf;
  plot(t,y(:,1));
end

function xp=F(t,x)
  xp=zeros(2,1); % since output must be a column vector
  xp(1)=x(2);
  xp(2)= -t*x(1) - exp(t)*x(2) + 3*sin(2*t);
end