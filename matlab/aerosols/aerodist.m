function aerodist
% setup the underlying  lognormal mass distribution
% see W\&H equation 6.8 for definitons
Ms=132.;
Mw=18.;
Rhol=1.e3;
Sigma=0.075;
rhoaero=1775.;
Rv=461.;
vanhoff=3;
Tinit=280.;

%distribution parameters
Mtot=.25e-9; %0.25 mug/m^3
logsd=1.7;
meandiam=0.2e-6;
meanmass=rhoaero*4./3.*pi.*(meandiam/2.).^3.
thelogmeanmass=log(meanmass)


%mass bins from 0.01 to 0.5 microns
%start with log10(mass)
drydiam= -2:.1:0.5;
%convert to meters
drydiam=10.^drydiam*1.e-6;
rhoaero=1775;
%convert aerosol volume to kg
themass=rhoaero*4./3.*pi.*(drydiam/2.).^3.;

%get the number distribtuion for themass vector
out=logDist(themass,Mtot,meanmass,logsd);
veclength=numel(out)

%now divide this by the bin mass, which gives
%the number of aerosols in each bin
thenum=out./themass;

%
% now work backwards to check our number distribution
%
dlogm= log(themass(2:end)) - log(themass(1:end-1));
disp('check total mass: ')
totmass=sum(out(2:end).*dlogm)
binnum=thenum(2:end).*dlogm;
disp('check total number: first')
totnum=sum(binnum)
thelogmeanmass=sum(out(2:end).*log10(themass(2:end)).*dlogm)/Mtot
integratedmeanmass=10^thelogmeanmass
fprintf('mean mass from N(m): %14.6e\n',calculatedmeanmass);

%plot the distribution
fh=figure(1)
semilogx(themass,out);
title('lognormal mass distribution');
xlabel('aerosol mass (kg)');
ylabel('m x mu(m) (kg/m^3/per log binwidth)');
print -depsc2 fig1.eps

%map masses to scrit
Scrit=mass2Scrit(themass);
Scrittrue=Scrit
themasstrue=log10(themass)

logS=log10(Scrit);
cumnum=fliplr(cumsum(fliplr(binnum)));


fh=figure(2);
size(Scrit)
size(cumsum(fliplr(binnum)))
Splot=Scrit(2:end);
loglog(Splot,cumnum,'r+');
%axis([0.01 2. 0 2.e8])
hold on
title('cumulative number distribution');
xlabel('Scrit');
ylabel('N (per m^3)');
disp('slope, intercept')
Svec=0.001:0.0001:.01;
Nvec=1.e9*Svec.^0.5;
drydiamvec=Scrit2drydiam(Svec);
massvec=Scrit2mass(Svec);
[rows Scount]=size(Nvec);
dN=(Nvec(2:Scount)-Nvec(1:Scount-1))
disp('Scrit, dry diam, massvec, Nvec')
[Svec;drydiamvec*1.e6;massvec*1.e19;Nvec*1.e-6]'
totmass=sum(massvec(2:Scount).*dN)*1.e9
totnum=-1*(Nvec(1) - Nvec(Scount))
plot(Svec,Nvec,'b-');
print -depsc2 fig2.eps
hold off;
%turn Svec into mass




a=(2.*Sigma)/(Rv*Tinit*Rhol);
bnomass=(vanhoff*Mw)/((4./3.)*pi*Rhol*Ms);
Scritnomass=((4*a^3.)./(27.*bnomass)).^0.5
[rows Scount]=size(Svec);
themass=zeros([rows Scount]);
for i = 1:Scount;
  themass(i)=(Scritnomass./Svec(i))^2.;
end
disp('here')
[log10(themass) ;log10(Svec)]'
[log10(themass) ; log10(Nvec)]'
dN=(Nvec(2:Scount)-Nvec(1:Scount-1))
disp('sum(dN)')
sum(dN)
disp('sum(dM)')
sum(dN.*themass(1:Scount-1))
a=(2.*Sigma)/(Rv*Tinit*Rhol);
b=(i*themass*Mw)/((4./3.)*pi*Rhol*Ms);
b=(i*Mw)/((4./3.)*pi*Rhol*Ms);
Scrit=((4*a^3.)./(27.*b)).^0.5
rcrit=(3*b/a).^0.5;

drydiam=[1.e-8 1.e-7];
Svec=drydiam2Scrit(drydiam);
Nvec=1.e9*Svec.^0.5
Nvec(1) - Nvec(2)
end


function out=logDist(themass,Amp,meanval,logsd);
   %out=logDist(themass,Amp,meanval,logsd)
   %return a vector with the mass (kg) in each bin specified by
   %themass, using the lognormal distribution given in the
   %aerosol notes.  Input: themass (kg), meanval(kg), logsd
   %(unitless)
   %output:  kg/logbinwidth, that is, out(i)*(log(themass(i+1)) - log(themass(i)))   
   %is the aerosol mass, in kg in bin(i), for aerosols with logmass
   %between log(themass(i)) and log(themass(i+1))
   theexp=exp(-0.5*(((log(themass) -log(meanval))/logsd)).^2.);
   out=Amp/(sqrt(2*pi)*logsd)*theexp;
end
   