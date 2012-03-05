 clear
 clc
 %plot a sounding from soundings.nc
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
%  fh=figure(1);
%  semilogy(temp,press)
%  hold on;
%  semilogy(dewpoint,press)
%  set(gca,'yscale','log','ydir','reverse');
%  ylim([400,1000]);
%  ylabel('press (hPa)')
%  xlabel('Temp (deg C)')
%  title('sounding 1')
%  hold off;
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
 

 pvec = 400:10:1000; %hPa
 for i=1:numel(pvec)
%   T_new=invtheta(thetaeVal, pvec(i)*100);
    [T_new, wv,wl] = tinvert_thetae(thetaeVal, 1, pvec(i)*100);
    xcoord_eq(i)=convertTempToSkew(T_new - c.Tc,pvec(i),skew);
 end
%  fprintf('ready to draw moist adiabat, thetae=%8.2f\n',thetaeVal);
semilogy(xcoord_eq,pvec,'k-','linewidth',3);
 ylim([400,1000.])
 hold off;
 
 %It does not matter if it is larger than 750 hPa
 %or smaller than 400 hPa because I will interpolate it
 ind = 350<press & press<800;
 newpress = press(ind);
 newtemp = temp(ind);
 
 press_int = 750:-1:400;
 temp_env=interp1(newpress, newtemp, press_int);
 
  for i=1:numel(press_int)
    [T_new, wv,wl] = tinvert_thetae(thetaeVal, 1, press_int(i)*100);
    temp_moist(i) = T_new-c.Tc;
  end
 
 temp_diff = temp_moist - temp_env;
 
 indPOS = temp_diff>=0;
 indNEG = temp_diff<0;
 
 temp_diff_cape = temp_diff(indPOS);
 newheight=press_int(indPOS);
 EL=newheight(1);
 temp_diff_cin = temp_diff(indNEG);

%Used the same method as exercise 8.1 to solve for CAPE and CIN
CAPE = c.Rd*mean(temp_diff_cape)*log(EL/400);
fprintf('thetae=%8.2f K\n',thetaeVal);
fprintf('Between 400-%3d hPa, CAPE=%8.2f J/kg\n', EL, CAPE);

if (EL == press_int(1))
    CIN = 0;
    fprintf('No Range, CIN=%8.2f J/kg\n', CIN);
else
    CIN = abs(c.Rd*mean(temp_diff_cin)*log(750/EL));
    fprintf('Between %3d-750 hPa, CIN=%8.2f J/kg\n', EL, CIN);
end


% $$$         double Mar-17-2011-00Z(dim_138, var_cols) ;
% $$$         double Mar-17-2011-12Z(dim_139, var_cols) ;
% $$$         double Mar-18-2011-00Z(dim_128, var_cols) ;
% $$$         double Mar-18-2011-12Z(dim_142, var_cols) ;
% $$$         double Mar-19-2011-00Z(dim_39, var_cols) ;
        