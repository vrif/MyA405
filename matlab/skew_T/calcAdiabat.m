function [pressVals,tempVals]=calcAdiabat(press0,thetae0,topPress)
    %[pressVals,tempVals]=calcAdiabat(press0,thetae0,topPress)
    %input: press0 (Pa),thetae0 (K), topPress (Pa)
    %output: pressVals (Pa), tempVals (K)
    c=constants;
    pressVals=linspace(press0,topPress,50);
    for i= 1:numel(pressVals)
      tempVals(i)=findTmoist(thetae0,pressVals(i));
      %fprintf('looping pressure: %8.2f\n',pressVals(i)*1.e-2);
    end
end
