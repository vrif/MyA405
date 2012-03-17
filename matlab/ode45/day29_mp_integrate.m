%find the rain rate for a dropsize distribution
%specified by a marshall-palmer distribution of 15 mm/hour
%first see how 
g=9.8;  %m/s^2
rho=1;  %kg/m^3 air density
rhol=1000.; %kg/m^3 liquid water density
N0=0.08*1.e6*1.e-1;%m^{-3} mm^{-1}
R=15; %mm/hr
dbar=1/(41.*R^(-2/9.)); %cm
dbar=dbar*10.; % mm
spacing=0.01; %mm
d=0:spacing:8;%mm
number=N0*exp(-d/dbar);% m^{-3} mm^{-1}
d=d*1.e-3; %meters
U=sqrt(rhol/rho*g*d); %m/s
R=sum(number*pi.*(d.^3)/6.*U*spacing); %m/s
R=R*1000*3600.; %mm/hour
out_mesg={'\nfor R=15 mm/hour integration gives %8.2f mm/hour\n'};
fprintf(strcat(out_mesg{:}),R);
