function [X,Y,Z,cdtr,V,sigma02,Qxx] = interf_calc_LS_code_GPS(PosSat,Dobs,X0)
%% function [X,Y,Z,dtr,V,sigma02,Qxx] = interf_calc_LS_code_GPS(PosSat,Dobs,X0)
%% Interface of function calc_LS_code
%% Only GPS, one pos and dte per epoch
%%
%% Clement Fontaine 2013-12-10
%%
%% Input : 
%% - PosSat : matrix containing satellite position [X,Y,Z] (m)
%% - Dobs : vector containing observations (m)
%% - X0 : initial values [X,Y,Z,dtr]
%%
%% Output : 
%% - X       : X (m) 
%% - Y       : Y (m) -> position in WGS84
%% - Z       : Z (m)
%% - cdtr    : cdtr (m)
%% - V       : residuals
%% - sigma02 : sigma^2 of compensation
%% - Qxx     : Var_covar matrix
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% constant
c = 299792458.0;


t = ones(size(Dobs)); % one date
ElevSat = (pi/2) * ones(size(Dobs)); % sigma = 2 m in calc_LS_code

sat_index = cell(size(Dobs,1),1);
for i = 1:size(Dobs,1)
	sat_index{i,1} = 'G';
end

if length(X0)<4
	X0 = [X0(:);zeros(4-length(X0),1)];
end

[result] = calc_LS_code(t,PosSat,Dobs,ElevSat,sat_index,X0);

% output
X = result.X;
Y = result.Y;
Z = result.Z;
cdtr = result.cdtr;
V = result.V;
Qxx = result.Qxx;
sigma02 = result.sigma02;


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
