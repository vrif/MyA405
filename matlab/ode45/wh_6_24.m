function wh_6_24
  tspan=0:0.1:2500;
  %start at cloudbase, with a 1 mm drop
  yinit=[0,1.e-3];
  [t,y]=ode45(@F,tspan,yinit);
  t=t/60.; %convert time to minutes
  height=y(:,1);
  radius=y(:,2)*1.e3;  %convert radius to mm
  figure(1);
  clf;
  plot(t,height);
  xlabel('time (minutes)');
  ylabel('height below cloud (m)');
  set(gca,'ydir','reverse');
  title('WH 6.24: drop position vs. time')
  grid on;
  figure(2);
  clf;
  plot(radius,height);
  set(gca,'ydir','reverse');
  xlabel('radius (mm)');
  ylabel('height below cloud (m)');
  title('WH 6.24: drop radius vs. height below cloudbase')
  set(gca,'ydir','reverse');
  grid on;
  %find the first index 5km below cloud base
  index_5k=find(abs(height - 5000.) < 0.5);
  index_5k=index_5k(1);
  therad=radius(index_5k);
  theheight=height(index_5k);
  thetime=t(index_5k);
  out_mesg={'\ndroplet is %8.2f meters below cloud base\n',...
          'at %8.3f mintues\n',...
          'with a radius of %8.3f mm\n\n'};
  fprintf(strcat(out_mesg{:}),theheight,thetime,therad);
end

function yp=F(t,y)
  yp=zeros(2,1); % since output must be a column vector
  %distance below cloud depends on radius
  yp(1)=6.e3*y(2);
  %radius depends on radius and saturation
  S=-0.4; %unitless
  Gl=700.e-12; %m^2/s  
  yp(2)= Gl*S/y(2);
end
