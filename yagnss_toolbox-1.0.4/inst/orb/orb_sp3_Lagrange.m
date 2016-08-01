function [X,Y,Z,dte,debug]=orb_sp3_Lagrange(sp3_data,mjd,constellation,PRN,degree)
%% function [X,Y,Z,dte,debug]=orb_sp3_Lagrange(sp3_data,mjd,constellation,PRN,degree)
%% Position and satellite clock error interpolation with a Lagrange polynomial
%%
%% Jacques Beilin - ENSG/DPTS - 2010-02-04
%% Clement Fontaine - 2013-11-12
%%
%% Input :
%% - sp3_data : structure created with function load_sp3.m
%% - mjd : Modified Julian Date of interpolation time
%% - constellation : 'G' for GPS, 'R' for Glonass and 'E' for Galileo
%% - PRN : satellite id
%% - degree : Lagrange polynomial degree
%%  
%% Output :
%% - X, Y, Z : position at the given mjd (m)
%% - dte : dte at the given mjd (s)
%% - debug : structure containing matrix extracted from sp3_data
%%
%% Position and dte are set to 0 if orbit is not computed
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extrapolation : allowed = 1, forbidden = 0
% WARNING : extrapolation with Lagrange polynomial is not recommended
extrapolate = 0;

X=0;
Y=0;
Z=0;
dte=0;
debug = cell(0);

m2 = ( degree + 1 ) / 2;

if ~(strcmp(constellation,'G') || strcmp(constellation,'R') || strcmp(constellation,'E'))
	tool_print_info(sprintf('Constellation %s is not implemented : coordinates and dte set to 0',constellation),2);
	return
end

if(PRN<1 || PRN>32)
	tool_print_info(sprintf('Satellite %s%02d does not exists : coordinates and dte set to 0',constellation,PRN),2);
	return
end

[orb,nl] = get_sp3(sp3_data,constellation,PRN);

% check if enougth epochs are available to interpolate a position with a Lagrange polynom
if degree > nl
	tool_print_info(sprintf('Satellite %s%02d : not enougth epochs to interpolate a position : coordinates and dte set to 0',constellation,PRN),2);
	return
end

% check presence of mjd into data period
if (extrapolate == 0)
	if mjd < orb(1,1) - 1/86400
		tool_print_info(sprintf('Satellite %s%02d : mjd is before first epoch of sp3 file (%.8f < %.8f) : coordinates and dte set to 0',constellation,PRN,mjd,orb(1,1)),2); 
		return
	end
	if mjd > orb(nl,1)
		tool_print_info(sprintf('Satellite %s%02d : mjd is after last epoch of sp3 file (%.8f > %.8f) : coordinates and dte set to 0',constellation,PRN,mjd,orb(nl,1)),2);  
		return
	end
end

% seek epoch just before mjd

pos = find(abs(mjd-orb(:,1))<15/1440);
pos = pos(1);

% side effects
if pos <= m2 % near file beginning
	first_index = 1;
elseif (pos < nl - m2) % normal case : (m + 1)/2 values around mjd
	first_index = pos - m2 + 1;
else % near file end
	first_index = nl - degree; 
end

% load temporary matrix
A = orb(first_index:first_index+degree,:);
debug.sp3_extract = A;

% interpolation and output
[Xs Ys Zs clock]=inter_Lagrange(A,mjd);

X = 1e3*Xs;
Y = 1e3*Ys;
Z = 1e3*Zs;
dte = 1e-6*clock;

end

function [Xs Ys Zs clock]=inter_Lagrange(sp3_extract,mjd)
%% [Xs Ys Zs clock]=inter_Lagrange(sp3_extract,mjd) : Lagrange interpolation in a matrix of positions and satellite clock errors
%%
%% Jacques Beilin - ENSG/DPTS - 2010-02-04
%% Clement Fontaine - 2013-10-22
%%
%% Input :
%% - sp3_extract : matrix extracted from sp3.G/R/E for one satellite (lines : epochs (degree + 1 values), columns : mjd X Y Z clk_error)
%% - mjd         : Modified Julian Date of interpolation time
%%  
%% Output :
%% - X : satellite position and clock error
%%

% détermination du degre en fonction de la taille de la matrice fournie

[nl,nc]=size(sp3_extract);
degree= nl-1;
Xs = 0;
Ys = 0;
Zs = 0;
clock = 0;
for (j=0:degree)
	Lj(j+1)=1.0;
	for (k=0:degree)
		if k~=j
			Lj(j+1) = Lj(j+1) * (mjd - sp3_extract(k+1,1)) / (sp3_extract(j+1,1) - sp3_extract(k+1,1));
		end
	end
	Xs = Xs + Lj(j+1) * sp3_extract(j+1,2);
	Ys = Ys + Lj(j+1) * sp3_extract(j+1,3);
	Zs = Zs + Lj(j+1) * sp3_extract(j+1,4);
	clock = clock + Lj(j+1) * sp3_extract(j+1,5);
end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
