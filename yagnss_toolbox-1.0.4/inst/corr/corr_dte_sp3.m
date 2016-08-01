function [dte] = corr_dte_sp3(sp3_data,const,PRN,mjd,degree);
%% [dte] = corr_dte_sp3(sp3_data,const,PRN,mjd,degree);
%% Satellite clock error from sp3.
%%
%% Clement Fontaine - 2013-12-04
%%
%% Input ,
%% - sp3_data : structure containing precise orbits, obtained from load_sp3
%% - const : costellation ('G' = GPS, 'R' = Glonass, 'E' = Galileo)
%% - PRN : satellite PRN
%% - mjd : modified Julian date
%% - degree : Lagrange interpolation degree
%%
%% Output
%% - dte : Satellite clock error (s)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Xs,Ys,Zs,dte] = orb_sat(sp3_data,const,PRN,mjd,degree);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
