function [dt_relat] = corr_dtrelat_sp3(sp3_data,const,PRN,mjd,delta,degree);
%% [dt_relat] = corr_dtrelat_sp3(sp3_data,const,PRN,mjd,delta,degree);
%% Relativistic correction calculation from sp3 data.
%%
%% Clement Fontaine - 2013-12-04
%%
%% Input 
%% - sp3_data : structure containing precise orbits, obtained from load_sp3
%% - const : costellation ('G' = GPS, 'R' = Glonass, 'E' = Galileo)
%% - PRN : satellite PRN
%% - mjd : modified Julian date
%% - delta : delta for derivation (1e-3 is ok)
%% - degree : Lagrange interpolation degree
%%
%% Output
%% - dt_relat : relativistic correction (s)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = 299792458.0;

delta_t = 1e-3; % delta for derivation

[Xs1,Ys1,Zs1,dte1] = orb_sat(sp3_data,const,PRN,mjd - delta_t/86400,degree);
[Xs2,Ys2,Zs2,dte2] = orb_sat(sp3_data,const,PRN,mjd + delta_t/86400,degree);

vec_Xs1 = [Xs1,Ys1,Zs1];
vec_Xs2 = [Xs2,Ys2,Zs2];
vec_Vx = (vec_Xs2 - vec_Xs1) / 2 / delta_t;
vec_Xs0 = [Xs1,Ys1,Zs1];
dt_relat = -2 * vec_Xs0 * vec_Vx' / c / c;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
