%comment
function updraft
   %plot a sounding 
    filename='springfield.nc';
    file_struct=nc_info(filename)
    c=constants;
    %
    % grap the first sounding pressure and temperature
    %
    sound_var = file_struct.Dataset(3).Name
    press=nc_varget(filename,sound_var,[0,0],[Inf,1]);
    temp=nc_varget(filename,sound_var,[0,2],[Inf,1]);
    dewpoint=nc_varget(filename,sound_var,[0,3],[Inf,1]);
    fh=figure(1);
    semilogy(temp,press)
    hold on;
    semilogy(dewpoint,press)
    set(gca,'yscale','log','ydir','reverse');
    ylim([400,1000]);
    ylabel('press (hPa)')
    xlabel('Temp (deg C)')
    title('sounding 1')
    hold off;
    figHandle=figure(2)
    skew=30.
    [figHandle,outputws,handlews]=makeSkew(figHandle,skew);
    xtemp=convertTempToSkew(temp,press,skew);    
    xdew=convertTempToSkew(dewpoint,press,skew);    
    semilogy(xtemp,press,'g-','linewidth',5);
    semilogy(xdew,press,'b-','linewidth',5);
    [xTemp,thePress]=ginput(1);
    Tclick=convertSkewToTemp(xTemp,thePress,skew);    
    thetaeVal=thetaes(Tclick + c.Tc,thePress*100.);
    [pressVals,tempVals]=calcAdiabat(thePress*100.,thetaeVal,400.e2);
    xTemp=convertTempToSkew(tempVals - c.Tc,pressVals*1.e-2,skew);
    fprintf('ready to draw moist adiabat, thetae=%8.2f\n',thetaeVal);
    semilogy(xTemp,pressVals*1.e-2,'r-','linewidth',5);
    ylim([400,1000.])
    xleft=convertTempToSkew(-20,1.e3,skew);
    xright=convertTempToSkew(25.,1.e3,skew);
    xlim([xleft,xright]);
    %
    % interpolator fails if two pressure values
    % are the same -- nudge them
    %
    newPress=nudgepress(press);
    interpTenv=@(pVals) interp1(newPress,temp,pVals);
    interpTdEnv=@(pVals) interp1(newPress,dewpoint,pVals);
    trytemp=interpTenv(pressVals*1.e-2);
    xTemp=convertTempToSkew(trytemp,pressVals*1.e-2,skew);
    semilogy(xTemp,pressVals*1.e-2,'b.','markersize',5);
    hold off;
    pressLevs=linspace(400,950,100)*1.e2;
    %
    % start integrating from 950 upwards
    %
    pressLevs=fliplr(pressLevs);
    TvDiffHandle=@(pVals) calcTvDiff(pVals,thetaeVal,interpTenv,interpTdEnv);
    for i=1:numel(pressLevs)
       Tvdiff(i)=TvDiffHandle(pressLevs(i));
    end
    figure(3)
    clf;
    plot(Tvdiff,pressLevs*0.01,'k-');
    set(gca,'ydir','reverse');
    ylabel('pressure (hPa)');
    xlabel('Virtual temperature difference (K)');
    title('Tvdiff vs. pressure')
    cumCAPE= -c.Rd*cumsum(Tvdiff(2:end).*diff(log(pressLevs)));
    figure(4)
    clf;
    plot(cumCAPE,pressLevs(2:end)*0.01,'k-');
    title('cumulative CAPE (J/kg) vs. pressure (hPa)')
    set(gca,'ydir','reverse');
    figure(5)
    clf;
    %
    % equate kinetic and potential energy toget maximum
    % updraft speed
    %
    maxvel=sqrt(2.*cumCAPE);
    plot(maxvel,pressLevs(2:end)*0.01,'k-');
    title('maximum updraft (m/s) vs. pressure (hPa)');
    set(gca,'ydir','reverse');
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

    
% $$$         double Mar-17-2011-00Z(dim_138, var_cols) ;
% $$$         double Mar-17-2011-12Z(dim_139, var_cols) ;
% $$$         double Mar-18-2011-00Z(dim_128, var_cols) ;
% $$$         double Mar-18-2011-12Z(dim_142, var_cols) ;
% $$$         double Mar-19-2011-00Z(dim_39, var_cols) ;
        