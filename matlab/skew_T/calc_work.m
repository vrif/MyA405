function calc_work
   [alphaCAB,pressCAB,alphaBC,pressBC]=calc_curves();
   plot_work(alphaCAB,pressCAB,alphaBC,pressBC);
end

function [alphaCAB,pressCAB,alphaBC,pressBC]=calc_curves
    c=constants;
    wtA=14.e-3;
    pressA=900.e2;
    tempA=25 + c.Tc;

    TdA=findTdwv(wtA,pressA);
    thetaeA=thetaep(TdA,tempA,pressA);
    wtB=wtA;
    pressB=700.e2;
    TdB=findTdwv(wtB,pressB);
    thetaeB=thetaeA;
    [tempB,wvB,wlB]=tinvert_thetae(thetaeB, wtB, pressB);
    wtC=wtA;
    wvA=wtA;
    pressC=900.e2;
    TdC=findTdwv(wtC,pressC);
    tempC=tempB;
    thetaeC=thetaep(TdC,tempC,pressC);
    pressLevs=linspace(700,900,60)*100.;
    for i=1:numel(pressLevs)
        thePress=pressLevs(i);
        [temp,wv(i),wl(i)]=tinvert_thetae(thetaeA, wtA, thePress);
        lineAB(i)=temp;
        rho=thePress/(c.Rd*temp);
        rhoAB(i)=rho;
    end
    [temp,wvC,wlC]=tinvert_thetae(thetaeC,wtC,900.e2);
    tempCA=linspace(tempC,tempA,100);
    rhoCA=pressA./(c.Rd*tempCA);
    press900Vec=NaN(size(rhoCA));
    press900Vec(:)=pressA;
    pressBC=pressLevs;
    rhoBC=pressBC/(c.Rd*tempB);
    alphaCA=1./rhoCA;
    alphaAB=1./rhoAB;
    alphaBC=1./rhoBC;
    alphaCAB=cat(2,alphaCA,fliplr(alphaAB));
    pressCAB=cat(2,press900Vec,fliplr(pressLevs));
    tempCA=linspace(tempC,tempA,100);
    rhoCA=pressA./(c.Rd*tempCA);
    press900Vec=NaN(size(rhoCA));
    press900Vec(:)=pressA;
    rhoBC=pressLevs/(c.Rd*tempB);
    Qin=c.cpd*(tempA - tempC) + c.lv0*(wvA - wvC);
    deltaS=c.cpd*(log(thetaeB) - log(thetaeC));
    Qout=tempB*deltaS;
    Qoutb= c.lv0*(wvB - wvC) - tempB*c.Rd*(log(pressB) - ...
                                           log(pressC));
    out_mesg={'Qout %9.3g (J/kg) via 2nd law \n',...
          'Qout %9.3g (J/kg) via first law\n',...
           'Qin %9.3g (J/kg)\n','net heating Qin - Qout %9.3g (J/kg)\n'};
    fprintf(strcat(out_mesg{:}),Qout,Qoutb,Qin,abs(Qin) - abs(Qout));
end    
    
function plot_work(alphaCAB,pressCAB,alphaBC,pressBC)
    figure(1);
    clf;
    interpCABfun=@(alphavals) interp1(alphaCAB,pressCAB,alphavals);
    interpBCfun=@(alphavals) interp1(alphaBC,pressBC,alphavals);
    interpAlphaVals=linspace(0.9,1.2,100);
    interpCAB=interpCABfun(interpAlphaVals);
    interpBC=interpBCfun(interpAlphaVals);
    plot(alphaCAB,pressCAB*1.e-2);
    ylim([700,1000.]);
    hold on;
    plot(alphaBC,pressBC*1.e-2);
    plot(interpAlphaVals,interpCAB*1.e-2,'r*');
    plot(interpAlphaVals,interpBC*1.e-2,'g*');
    start=alphaBC(end);
    stop=alphaBC(1);
    workBC=quad(interpBCfun,start,stop);
    workCAB=quad(interpCABfun,start,stop);
    netWork=abs(workCAB) - abs(workBC);
    out_mesg={'workBC %9.3g (J/kg)\n',...
          'workCAB %9.3g (J/kg)\n',...
          'net work done by system %9.3g (J/kg)\n'};
    fprintf(strcat(out_mesg{:}),workBC,workCAB,abs(workCAB) - abs(workBC));
    title('volume - temperature plot');
    hold off;
end        

