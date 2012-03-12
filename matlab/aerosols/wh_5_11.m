function wh_5_11
%distribution parameters
Stot=1300.; %micron^2/cm^3
logsd=1.1;
meanDiam=0.1; %micron

%span the range of the first peak figure 5_11
drydiam= -3:.05:1.3;
drydiam=10.^drydiam;  %microns
%get the lognormal surface area distribtuion for the diameters
sdist=logDist(drydiam,Stot,meanDiam,logsd);
%integrate this to check area
sbin=drydiam(2:end).*sdist(2:end); %surface area in the logDiam bin
check=sum(sbin.*diff(log(drydiam)));
fprintf('check total mass is 1300 micron^2/cm^3: %8.2f\n',check)
%number in each log bin
%plot it
fh=figure(1);
clf;
semilogx(drydiam(2:end),sbin);
title('lognormal surface area distribution');
xlabel('aerosol diameter (microns)');
ylabel('D \times sdist (microns^2/cm^3/per log binwidth)');
print -dpng -r150 Dsdist.png;
print -depsc  Dsdist.eps;

%find the number distribution with units number/cm^3/log binwidth
thearea=4.*pi*(drydiam/2.).^2.;
nbin=sbin./thearea(2:end);
Ntot=sum(nbin.*diff(log(drydiam)));
fprintf('the total number/cm^3: %8.2f\n',Ntot)

fh=figure(2);
clf;
loglog(drydiam(2:end),nbin);
title('number distribution');
xlabel('aerosol diameter (microns)');
ylabel('D \times ndist (number/cm^3/per log binwidth)');
print -dpng -r150 Dndist.png;
print -depsc  Dndist.eps;



end

% $$$ 
% $$$ %plot number distribution
% $$$ fh=figure(2);
% $$$ clf;
% $$$ semilogx(themass,themass.*ndist);
% $$$ title('n(m) for lognormal mass distribution');
% $$$ xlabel('aerosol mass (kg)');
% $$$ ylabel('m \times ndist (number/m^3/per log binwidth)');
% $$$ print -dpng -r150 mndist.png;
% $$$ print -depsc mndist.eps;
% $$$ 
% $$$ ndiam=ndist(2:end).*diff(themass)./diff(drydiam);
% $$$ 
% $$$ fh=figure(3);
% $$$ clf;
% $$$ %note the units: microns for diameter, #/cc/logbinwidth for ndiam
% $$$ semilogx(drydiam(2:end).*1.e6,drydiam(2:end).*ndiam);
% $$$ title('ndiam(D) for lognormal mass distribution');
% $$$ xlabel('aerosol diameter (microns)');
% $$$ ylabel('D \times ndist (number/cm^3/per log binwidth)');
% $$$ print -dpng -r150 Dndiam.png;
% $$$ print -depsc Dndiam.eps;
% $$$ Ntot=sum(ndiam.*diff(drydiam));
% $$$ %1.e-6 converts from m^3 to cm^3
% $$$ fprintf('check the total number/cm^3: %8.2f\n',Ntot*1.e-6)
% $$$ 
% $$$ 
% $$$ Nbins=ndiam.*diff(drydiam);
% $$$ N=cumsum(Nbins);
% $$$ N=fliplr(N);
% $$$ fh=figure(4);
% $$$ clf;
% $$$ %note the units: microns for diameter, #/cc/logbinwidth for ndiam
% $$$ semilogx(drydiam(2:end).*1.e6,N);
% $$$ title('ndiam(D) for lognormal mass distribution');
% $$$ xlabel('aerosol diameter (microns)');
% $$$ ylabel('D \times ndist (number/cm^3/per log binwidth)');
% $$$ print -dpng -r150 N.png;
% $$$ print -depsc N.eps;
% $$$ Ntot=sum(N.*diff(drydiam));
% $$$ %1.e-6 converts from m^3 to cm^3
% $$$ fprintf('check the total number/cm^3: %8.2f\n',Ntot*1.e-6)
% $$$ 
% $$$ 
% $$$ end
% $$$ 
% $$$ 

function sdist=logDist(theDiam,Stot,meanDiam,logsd);
   %sdist=logDist(theDiam,Stot,meanval,logsd)
   %return a vector with the surface are (m^2/m^3) in each bin specified by
   %theDiam, using the lognormal distribution given in the
   %aerosol notes.  Input: theDiam (m), meanval(m), logsd
   %(unitless)
   %output:  microns^2/cm^3/micron binwidth, 
   %if Stot is in micons^2/cm^3 and theDiam is in microns
   %that is, theDiam(i)*sdist(i)*(log(theDiam(i+1)) - log(theDiam(i))) 
   % = sdist(i)*(theDiam(i+1) - theDiam(i))
   %is the aerosol surface area concentration, in m^2/m^3 in bin(i), for aerosols with diameters
   %between theDiam(i)) and theDiam(i+1)
   theexp=exp(-0.5*(((log(theDiam) -log(meanDiam))/logsd)).^2.);
   sdist=Stot./theDiam/(sqrt(2*pi)*logsd).*theexp;
end
   