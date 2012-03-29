function day31

   %plot a sounding 
    filename='littlerock.nc';
    fprintf('reading file: %s\n',filename);
    file_struct=nc_info(filename);
    c=constants;
    %
    % grab the March 2 12Z sounding
    %
    sound_var = file_struct.Dataset(4).Name;
    fprintf('found sounding: %s\n',sound_var);
    press=nc_varget(filename,sound_var,[0,0],[Inf,1]);
    temp=nc_varget(filename,sound_var,[0,2],[Inf,1]);
    dewpoint=nc_varget(filename,sound_var,[0,3],[Inf,1]);
    height=nc_varget(filename,sound_var,[0,1],[Inf,1]);
    newheight = nudgeHeight(height);
%     fh=figure(1);
%     semilogy(temp,press);
%     hold on;
%     semilogy(dewpoint,press);
%     set(gca,'yscale','log','ydir','reverse');
%     ylim([400,1000]);
%     ylabel('press (hPa)');
%     xlabel('Temp (deg C)');
%     title('sounding 1');
%     hold off;
%     figHandle=figure(2);
%     skew=30.;
%     [figHandle,outputws,handlews]=makeSkew(figHandle,skew);
%     xtemp=convertTempToSkew(temp,press,skew);
%     xdew=convertTempToSkew(dewpoint,press,skew);
%     semilogy(xtemp,press,'g-','linewidth',5);
%     semilogy(xdew,press,'b-','linewidth',5);
%     %use lowest sounding level for adiabat
    thetaeVal=thetaes(temp(1) + c.Tc,press(1)*100.);
    [pressVals,tempVals]=calcAdiabat(press(1)*100.,thetaeVal,400.e2);
%     xTemp=convertTempToSkew(tempVals - c.Tc,pressVals*1.e-2,skew);
%     semilogy(xTemp,pressVals*1.e-2,'r-','linewidth',5);
%     ylim([400,1000.]);
%     xleft=convertTempToSkew(-20,1.e3,skew);
%     xright=convertTempToSkew(25.,1.e3,skew);
%     xlim([xleft,xright]);
    %
    % interpolator fails if two pressure values
    % are the same -- nudge them
    %
    newPress=nudgepress(press);
    interpTenv=@(pVals) interp1(newPress,temp,pVals);
    interpTdEnv=@(pVals) interp1(newPress,dewpoint,pVals);
    trytemp=interpTenv(pressVals*1.e-2);
%     xTemp=convertTempToSkew(trytemp,pressVals*1.e-2,skew);
%     semilogy(xTemp,pressVals*1.e-2,'b.','markersize',5);
%     hold off;
    pressLevs=linspace(100,press(1),100)*1.e2;
    %
    % starting integrating from first sounding level
    %
    pressLevs=fliplr(pressLevs);
    TvDiffHandle=@(pVals) calcTvDiff(pVals,thetaeVal,interpTenv,interpTdEnv);
    for i=1:numel(pressLevs)
       Tvdiff(i)=TvDiffHandle(pressLevs(i));
    end
    
    TvDiffHandle=@(pVals) calcTvDiff(pVals,thetaeVal,interpTenv,interpTdEnv);
    pressLevs=linspace(100,press(1),100)*1.e2;
    %
    % starting Buoyancy from first sounding
    %
    pressLevs=fliplr(pressLevs);
    BHandle=@(pVals) calcBuoyancy(pVals,thetaeVal,interpTenv,interpTdEnv);
    for i=1:numel(pressLevs)
       B(i)=BHandle(pressLevs(i));
    end
    
%     figure(2);
%     clf;
%     plot(Tvdiff,pressLevs*0.01,'k-');
%     set(gca,'ydir','reverse');
%     ylabel('pressure (hPa)');
%     xlabel('Virtual temperature difference (K)');
%     title('Tvdiff vs. pressure');
    figure(1);
    clf;
    plot(B,pressLevs*0.01,'k-');
    set(gca,'ydir','reverse');
    ylabel('pressure (hPa)');
    xlabel('Buoyancy (m/s^2)');
    title('Buoyancy vs Pressure');
    
    yinit=[0,0];
    derivs=@(t,y) F(t,y,thetaeVal,interpTenv,interpTdEnv);
    [t,y]=ode45(derivs,pressLevs,yinit);
    
    figure(2);
    clf;
    plot(y(:,2),pressLevs*0.01,'k-');
    set(gca,'ydir','reverse');
    ylabel('pressure (hPa)');
    xlabel('w vertical velocity (m/s)');
    title('W, Vertical Velocity vs. Pressure');
    
    pressLevs=linspace(200,press(1),100)*1.e2;
    %
    % starting integrating from first sounding level
    %
    pressLevs=fliplr(pressLevs);
    TvDiffHandle=@(pVals) calcTvDiff(pVals,thetaeVal,interpTenv,interpTdEnv);
    for i=1:numel(pressLevs)
       Tvdiff(i)=TvDiffHandle(pressLevs(i));
    end
    cumCAPE= -c.Rd*cumsum(Tvdiff(2:end).*diff(log(pressLevs)));
    figure(4);
    clf;
    plot( cumCAPE,pressLevs(2:end)*0.01,'k-');
    title('cumulative CAPE (J/kg) vs. pressure (hPa)');
    set(gca,'ydir','reverse'); 
    
    figure(5);
    maxvel=sqrt(2.*cumCAPE);
    plot(maxvel,pressLevs(2:end)*0.01,'k-');
    title('maximum updraft (m/s) vs. pressure (hPa)');
    set(gca,'ydir','reverse');
    derivs;
  
    
end

function yp=F(t,y,theta0,interpTenv,interpTdEnv)
  yp=zeros(2,1); % since output must be a column vector
  BHandle=@(pVals) calcBuoyancy(pVals,theta0,interpTenv,interpTdEnv);
  %first variable is w updraft in m/s
  yp(2)=BHandle(t);
  %second variable is height in m
  yp(1)= yp(2);
end

function newPress=nudgepress(pressVec)
    %if two balloon pressure levels are idential
    %add a factor of 0.1% to the second one
    %so interpolation will work
    newPress=pressVec;
    hit=find(abs(diff(newPress)) < 1.e-8);
    newPress(hit+1)=pressVec(hit) + 1.e-3*pressVec(hit);
end

function TvDiff=calcTvDiff(press,thetae0,interpTenv,interpTdEnv)
    %calcTvDiff(press,thetae0,interpTenv,interpTdenv)
    %input: press (Pa), thetae0 (K), plus function handles for T,Td soundings
    %output: TvDiff (K)
    %neglect liquid water loading in the virtual temperature
    c=constants;
    Tcloud=findTmoist(thetae0,press);
    wvcloud=wsat(Tcloud,press);
    Tvcloud=Tcloud*(1. + c.eps*wvcloud);
    Tenv=interpTenv(press*1.e-2) + c.Tc;
    Tdenv=interpTdEnv(press*1.e-2) + c.Tc;
    wvenv=wsat(Tdenv,press);
    Tvenv=Tenv*(1. + c.eps*wvenv);
    TvDiff=Tvcloud - Tvenv;
end


function B=calcBuoyancy(press,thetae0,interpTenv,interpTdEnv)
    TvDiffHandle=@(pVals) calcTvDiff(pVals,thetae0,interpTenv,interpTdEnv);
    Tv_hat = TvDiffHandle(press);
    
    c=constants;
    Tenv=interpTenv(press*1.e-2) + c.Tc;
    Tdenv=interpTdEnv(press*1.e-2) + c.Tc;
    wvenv=wsat(Tdenv,press);
    Tvenv=Tenv*(1. + c.eps*wvenv);
    Tv_bar = Tvenv;
    
    B = 9.8*(Tv_hat/Tv_bar);
end

function newHeight=nudgeHeight(zVec)
   %if two balloon height levels are idential or descend
   %add a factor of 0.1% to the second one
   %so interpolation will work
   newHeight=zVec;
   hit=find(diff(newHeight) <= 0);
   newHeight(hit+1)=zVec(hit) + 1.e-3*zVec(hit);
end


