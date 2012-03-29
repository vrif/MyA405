function littlerock_stats
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
    direct=nc_varget(filename,sound_var,[0,6],[Inf,1]);
    speed=nc_varget(filename,sound_var,[0,7],[Inf,1]);
    fh=figure(1);
    semilogy(temp,press);
    hold on;
    semilogy(dewpoint,press);
    set(gca,'yscale','log','ydir','reverse');
    ylim([200,1000]);
    ylabel('press (hPa)');
    xlabel('Temp (deg C)');
    title('sounding 1');
    hold off;
    figHandle=figure(2);
    skew=30.;
    [figHandle,outputws,handlews]=makeSkew(figHandle,skew);
    xtemp=convertTempToSkew(temp,press,skew);
    xdew=convertTempToSkew(dewpoint,press,skew);
    semilogy(xtemp,press,'g-','linewidth',5);
    semilogy(xdew,press,'b-','linewidth',5);
    %use 900 hPa sounding level for adiabat
    p900_level=find(abs(900 - press) < 2.);
    thetaeVal=thetaep(dewpoint(p900_level) + c.Tc,temp(p900_level) + c.Tc,press(p900_level)*100.);
    [pressVals,tempVals]=calcAdiabat(press(p900_level)*100.,thetaeVal,200.e2);
    xTemp=convertTempToSkew(tempVals - c.Tc,pressVals*1.e-2,skew);
    semilogy(xTemp,pressVals*1.e-2,'r-','linewidth',5);
    ylim([200,1000.]);
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
    interpDirec=@(pVals) interp1(newPress,direct,pVals);
    interpSpeed=@(pVals) interp1(newPress,speed,pVals);
    trytemp=interpTenv(pressVals*1.e-2);
    xTemp=convertTempToSkew(trytemp,pressVals*1.e-2,skew);
    semilogy(xTemp,pressVals*1.e-2,'b.','markersize',5);
    hold off;
    figure(3);
    clf;
    pressLevs=linspace(200,press(1),100)*1.e2;
    %
    % starting integrating from first sounding level
    %
    pressLevs=fliplr(pressLevs);
    TvDiffHandle=@(pVals) calcTvDiff(pVals,thetaeVal,interpTenv,interpTdEnv);
    for i=1:numel(pressLevs)
       Tvdiff(i)=TvDiffHandle(pressLevs(i));
    end
    plot(Tvdiff,pressLevs*0.01,'k-');
    set(gca,'ydir','reverse');
    ylabel('pressure (hPa)');
    xlabel('Virtual temperature difference (K)');
    title('Tvdiff vs. pressure');
    cumCAPE= -c.Rd*cumsum(Tvdiff(2:end).*diff(log(pressLevs)));
    figure(4);
    clf;
    plot(cumCAPE,pressLevs(2:end)*0.01,'k-');
    title('cumulative CAPE (J/kg) vs. pressure (hPa)');
    set(gca,'ydir','reverse');
    figure(5);
    clf;
    %
    % equate kinetic and potential energy toget maximum
    % updraft speed
    %
    maxvel=sqrt(2.*cumCAPE);
    plot(maxvel,pressLevs(2:end)*0.01,'k-');
    title('maximum updraft (m/s) vs. pressure (hPa)');
    set(gca,'ydir','reverse');
%
% find storm indices
%
%  lifted index 
    thetaeVal=thetaep(dewpoint(1) + c.Tc,temp(1) + c.Tc,press(1)*100.);
    wT=wsat(dewpoint(1) + c.Tc,press(1)*100.);
    Tadia_500=tinvert_thetae(thetaeVal,wT,500.e2);
    Temp_500=interpTenv(500.) + c.Tc;
    lifted_index= Temp_500 - Tadia_500;
% total totals = vertical totals plus cross totals
    Temp_850=interpTenv(850.) + c.Tc;
    dew_850 = interpTdEnv(850.) + c.Tc;
    TT_index=Temp_850 + dew_850 - 2*Temp_500;
%  Sholwater
    thetaeVal=thetaep(dew_850,Temp_850,850.*100.);
    wT=wsat(dew_850,850*100.);
    Tadia_500=tinvert_thetae(thetaeVal,wT,500.e2);
    sholwater=Temp_500 - Tadia_500;
%  SWEAT
    speed_850=interpSpeed(850.);
    speed_500=interpSpeed(500.);
    dir_850=interpDirec(850.);
    dir_500=interpDirec(500.);
    angle=(dir_500 - dir_850);
    shear=125.*(sind(angle) + 0.2);
    sweat=12.*(dew_850 - c.Tc) + 20.*(TT_index - 49.) +  ...
           2.*speed_850 + speed_500 + shear;
    out_mesg={'\nCAPE is %10.4f (J/kg)\n', ...
              '\nLifted index= %10.4f (K)\n', ...
              '\nTotal Totals=%10.4f (K)\n', ...
              '\nSholwater=%10.4f (K)\n', ...
              '\nSWEAT =%10.4f\n'};
    fprintf(strcat(out_mesg{:}),cumCAPE(end),...
            lifted_index,TT_index,sholwater,sweat);
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
