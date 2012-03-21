c=constants();
pressA=1.e5;
tempA=15 + c.Tc;
pressC=0.7e5;
tempC=5 + c.Tc;
thetaA=theta(tempA,pressA);
thetaC=theta(tempC,pressC);
thetaB=thetaC;
tempB=tempA;
term1=(thetaB/tempB)^(c.cpd/c.Rd);
term1=1./term1;
pressB=term1*1.e5;
tempD=tempC;
thetaD=thetaA;
term1=(thetaD/tempD)^(c.cpd/c.Rd);
term1=1./term1;
pressD=term1*1.e5;
figHandle=figure(1);
clf;
skew=30;
[figureHandle,outputws,handlews]=makeSkew(figHandle,skew);
xtempA=convertTempToSkew(tempA - c.Tc,pressA*0.01,skew);
xtempB=convertTempToSkew(tempB - c.Tc,pressB*0.01,skew);
xtempC=convertTempToSkew(tempC - c.Tc,pressC*0.01,skew);
xtempD=convertTempToSkew(tempD - c.Tc,pressD*0.01,skew);
%plot(xtempA,pressA*0.01,'bd','markerSize',12,'markerFaceColor','b');
text(xtempA,pressA*0.01,'A','edgecolor','b','fontweight','bold','fontsize',22,'color','b');
text(xtempB,pressB*0.01,'B','edgecolor','b','fontweight','bold','fontsize',22,'color','b');
text(xtempC,pressC*0.01,'C','edgecolor','b','fontweight','bold','fontsize',22,'color','b');
text(xtempD,pressD*0.01,'D','edgecolor','b','fontweight','bold','fontsize',22,'color','b');
print -dpdf forward_carnot.pdf
print -dpng -r200 forward_carnot.png
thetaB=thetaA;
tempB=tempC;
term1=(thetaB/tempB)^(c.cpd/c.Rd);
term1=1./term1;
pressB=term1*1.e5;
tempD=tempA;
thetaD=thetaC;
term1=(thetaD/tempD)^(c.cpd/c.Rd);
term1=1./term1;
pressD=term1*1.e5;
figHandle=figure(2);
clf;
skew=30;
[figureHandle,outputws,handlews]=makeSkew(figHandle,skew);
xtempA=convertTempToSkew(tempA - c.Tc,pressA*0.01,skew);
xtempB=convertTempToSkew(tempB - c.Tc,pressB*0.01,skew);
xtempC=convertTempToSkew(tempC - c.Tc,pressC*0.01,skew);
xtempD=convertTempToSkew(tempD - c.Tc,pressD*0.01,skew);
%plot(xtempA,pressA*0.01,'bd','markerSize',12,'markerFaceColor','b');
text(xtempA,pressA*0.01,'A','edgecolor','b','fontweight','bold','fontsize',22,'color','b');
text(xtempB,pressB*0.01,'B','edgecolor','b','fontweight','bold','fontsize',22,'color','b');
text(xtempC,pressC*0.01,'C','edgecolor','b','fontweight','bold','fontsize',22,'color','b');
text(xtempD,pressD*0.01,'D','edgecolor','b','fontweight','bold','fontsize',22,'color','b');
print -dpdf backward_carnot.pdf
print -dpng -r200 backward_carnot.png
print -depsc  backward_carnot_clip.eps

%Calculate the heat input and output

%exact change in entropy
deltaS=c.cpd*(log(thetaC) - log(thetaB));
%approx change in entropy
deltaSapprox=c.cpd*(thetaC - thetaB)/thetaC;
%heat input inside the house
qin=tempB*deltaS;
%heat exausted outside the house
qout=tempA*deltaS;
work_done=qin - qout;


out_mesg={'\nHeat absorbed during expansion (qin) %8.3f (J/kg)%%\n',...
          'Heat expelled during compression (qout) %8.3f (J/kg)\n',...
           'work required to complete the cycle (work_done) %8.3f (J/kg)\n',...
           '(this is also the net energy removed from the room)\n'};
fprintf(strcat(out_mesg{:}),qin,qout,work_done);

%check this using enthalpy per the carnot_refrigerator notes

qin_h= c.cpd*(tempC - tempB) - c.Rd*tempB*(log(pressC) - log(pressB));
qout_h= c.cpd*(tempD - tempA) - c.Rd*tempD*(log(pressD) - log(pressA));

out_mesg={'\n\nEnthalpy check: heat absorbed during expansion (qin) %8.3f (J/kg)%%\n',...
          'Enthalpy check: Heat expelled during compression (qout) %8.3f (J/kg)\n'};

fprintf(strcat(out_mesg{:}),qin_h,qout_h);
