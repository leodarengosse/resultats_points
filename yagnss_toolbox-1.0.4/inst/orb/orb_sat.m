function [X,Y,Z,dte,debug] = orb_sat(data,const,PRN,mjd,degree);
%% function [X,Y,Z,dte,debug] = orb_sat(data,const,PRN,mjd,degree);
%%
%% Calculates ECEF GPS, Galileo and Glonass satellite coordinates and 
%% satellite clock error from ephemeris structure or sp3_data
%%
%% Jacques Beilin - ENSG/DPTS - 2010-02-25
%% Clement Fontaine - 2013-11-12
%%
%% Input :
%% - data : Eph structure obtained from get_ephemeris or sp3_data obtained 
%%   from load_sp3 function
%% - const : constellation id ('G' for GPS, 'R' for Glonass and 'E' for Galileo) 
%% - PRN : satellite id in constellation
%% - mjd : modified Julian day
%% - degree : degree for Lagrange interpolation (optional)
%%
%% WARNING : if data is a Eph structure, constellation and PRN used come
%%           from Eph structure (Eph.const and Eph.PRN)
%%
%% Output 
%% - X,Y,Z : cartesian coordinates : 
%%      - if Eph structure : 
%%			- WGS84 for GPS 
%%			- GTRF for Galileo 
%%			- WGS84 for Glonass (PZ90 to WGS84 in pos_Glonass)
%%      - else : sp3 coordinate system (IGS08)
%% - dte : satellite clock offset
%% - debug : debug structure with all results.
%%
%% If orbit is not computed, position and dte are set to 0
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X=0;
Y=0;
Z=0;
dte=0;
debug = cell(0);


% Select case : sp3 or brdc
if(isfield(data,'Gnb') || isfield(data,'Rnb') || isfield(data,'Enb')) % sp3_data structure + test if structure is empty

	sp3_data = data;
	flag = 1;

elseif(isfield(data,'mjd')) % Eph structure + test if structure is empty

	Eph = data;
	flag = 2;

else
	tool_print_info('Eph structure is empty : coordinates and dte set to 0',2);
	flag = 0; % structure id empty
	return;

end


% Case 1 : sp3 file

if flag == 1

	% Lagrange interpolation
	[X,Y,Z,dte,debug]=orb_sp3_Lagrange(sp3_data,mjd,const,PRN,degree);
	
	% Center of mass 2 center of phase
	% TODO
	

% case 2 : Eph file

elseif flag == 2

	% test : GPS and Galileo (keplerian elements) or Glonass (numerical integration)

	if (strcmp(Eph.const,'G') || strcmp(Eph.const,'E')) % GPS and Galileo
	
		[X,Y,Z,dte,debug] = orb_from_eph(Eph,mjd);
		
	elseif strcmp(Eph.const,'R'); % Glonass
	
		[t,X,Y,Z,VXs,VYs,VZs,dte,debug] = orb_from_RK4(Eph,mjd);

	else
		tool_print_info('Constellation not implemented',2);
		return
		
	end

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
