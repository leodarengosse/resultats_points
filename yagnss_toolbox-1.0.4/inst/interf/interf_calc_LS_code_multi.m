function [X,Y,Z,cdtr,cGGTO,cGPGL,V,sigma02,Qxx] = interf_calc_LS_code_multi(PosSat,Dobs,sat_index,X0)
%% function [X,Y,Z,cdtr,cGGTO,cGPGL,V,sigma02,Qxx] = interf_calc_LS_code_multi(PosSat,Dobs,sat_index,X0)

%% Interface of function calc_LS_code
%% Only GPS, one pos and dte per epoch
%%
%% Clement Fontaine 2013-12-10
%%
%% Input : 
%% - PosSat : matrix containing satellite position [X,Y,Z] (m)
%% - Dobs : vector containing observations (m)
%% - sat_index : vector containing constellation type (1 = GPS, 2 = Galileo, 3 = Glonass)
%% - X0 : initial values [X,Y,Z,dtr,cGGTO,cGPGL]
%%
%% Output : 
%% - X       : X (m) 
%% - Y       : Y (m) -> position in WGS84
%% - Z       : Z (m)
%% - cdtr    : c * dtr (m)
%% - cGGTO   : c * GPS to Galileo Time Offset (m)
%% - cGPGL   : c * Glonass to GPSTime (m)
%% - V       : residuals
%% - sigma02 : sigma^2 of compensation
%% - Qxx     : Var_covar matrix
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% constant
c = 299792458.0;

if length(X0)<6
	X0 = [X0(:);zeros(6-length(X0),1)];
end

t = ones(size(Dobs)); % one date
ElevSat = (pi/2) * ones(size(Dobs)); % sigma = 2 m in calc_LS_code
sat_index_temp = cell(size(Dobs,1),2);
for i = 1:size(Dobs,1)
	if (sat_index(i) == 1)
		sat_index_temp{i,1} = 'G';
	elseif (sat_index(i) == 2)
		sat_index_temp{i,1} = 'E';
	elseif (sat_index(i) == 3)
		sat_index_temp{i,1} = 'R';
	end
end

sat_index = sat_index_temp;

[result] = calc_LS_code(t,PosSat,Dobs,ElevSat,sat_index,X0);

% output
X = result.X;
Y = result.Y;
Z = result.Z;
cdtr = result.cdtr;
cGPGL = result.cGPGL;
cGGTO = result.cGGTO;
V = result.V;
Qxx = result.Qxx;
sigma02 = result.sigma02;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
