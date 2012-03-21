function TvDiff=calcTvDiff(press,thetae0,interpTenv,interpTdEnv)
    %calcTvDiff(press,thetae0,interpTenv,interpTdenv)
    %input: press (Pa), thetae0 (K), plus function handles for T,Td soundings
    %output: TvDiff (K)
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
