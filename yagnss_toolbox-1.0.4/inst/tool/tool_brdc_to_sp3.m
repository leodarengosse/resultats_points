function [sp3_header,sp3_data] = tool_brdc_to_sp3(NAV_header,NAV_data,mjd_min,mjd_max,step)
%% function [sp3_header,sp3_data] = tool_brdc_to_sp3(NAV_header,NAV_data,mjd_min,mjd_max,step)
%%
%% Compute orbit of satellites which are present in NAV_data, in order to set a sp3 structure
%%
%% Clement Fontaine - 2013-11-13
%%
%% Input : 
%% - NAV_header, NAV_data : navigation message (from load_rinex_n)
%% - mjd_min, mjd_max : min and max of sp3 (mjd)
%% - step : step in seconds
%%
%% Output : 
%% - sp3_header, sp3_data : sp3 formated orbits and dtr (see load_sp3 for details)
%%
%%   Empty structure is returned if NAV_header format isn't good
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Output
sp3_data = cell(0);
sp3_header = cell(0);

% Test NAV structure
if ~isfield(NAV_header,'GPSA')
	tool_print_info('ERROR : input a NAV structure (see function load_rinex_n)',3)
	return;
end

% Nuber of epochs
N_epoch = length(mjd_min:step/86400:mjd_max);

gpst = mjd_t(mjd_min);

% SP3 header
sp3_header.Version = NAV_header.VERSION;
sp3_header.Flag = '';
sp3_header.Date = sprintf('%4d %02d %02d %02d %02d %f', gpst.yyyy, gpst.mon, gpst.dd, gpst.hh, gpst.min, gpst.sec);
sp3_header.Number_of_Epochs= N_epoch;

sp3_header.Data_Used = 'EPH';
sp3_header.Coordinate_Sys = 'WGS84 for GPS and GLO, GTRF for GAL';
sp3_header.Orbit_Type = NAV_header.TYPE;
sp3_header.Agency = '';
%~ sp3_header.Leap_seconds = NAV_header.LEAP_SECONDS;

sp3_header.wk = gpst.wk;
sp3_header.sow = gpst.wsec; 
sp3_header.Epoch_Interval = step; 
sp3_header.mjd = gpst.mjd;
sp3_header.Fractional_Day = ''; 

% SP3_data
const = 'GRE'; % GPS, GLONASS and GALILEO supported

for i_const = 1:length(const) % constellation

	sp3_data_const = zeros(32,5,N_epoch);
	nb_sat = zeros(32,1);

	for PRN = 1:32 % PRN

		ind = 1;
		
		for i = mjd_min:step/86400:mjd_max
				
			[Eph]=get_ephemeris(NAV_header,NAV_data,const(i_const),PRN,i);
			[Xs,Ys,Zs,dte,debug] = orb_sat(Eph,const(i_const),PRN,i);
			sp3_data_const(PRN,:,ind) = [i,Xs,Ys,Zs,dte];
			
			if (Xs~=0) % test if Xs exists (= non 0)
				nb_sat(PRN) = nb_sat(PRN) + 1;
			end
			
			ind = ind+1;
		
		end
		
	end
	
	sp3_data.(const(i_const)) = sp3_data_const;
	sp3_data.(strcat(const(i_const),'nb')) = nb_sat;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end
