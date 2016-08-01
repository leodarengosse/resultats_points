function [Obs_base, Obs_rover] = get_obs_in_common(RNX_header_base,RNX_data_base, RNX_header_rover, RNX_data_rover, constellation, PRN, epoch)
%% function [Obs_base, Obs_rover] = get_obs_in_common(RNX_header_base,RNX_data_base, RNX_header_rover, RNX_data_rover, constellation, PRN, epoch)
%%
%% Get obs data in common between 2 stations, for one satellite at a given epoch (GPS, Glonass and Galileo supported)
%%
%% Clement Fontaine 2014-01-27
%%
%% Input :
%% - RNX_header_base : structure containing Rinex header for base station
%% - RNX_data_base : matrix containing data for base station
%% - RNX_header_rover : structure containing Rinex header for rover station
%% - RNX_data_rover : matrix containing data for rover station
%%   RNX_header and RNX_data are set up with function load_rinex_o
%% - constellation : 'G' = GPS, 'R' = GLONASS, 'E' = Galileo
%% - PRN : satellite id
%% - epoch
%%
%% Output : 
%% 	- Obs_base, Obs_rover : structures containing informations
%%    - constellation
%%    - PRN
%%    - epoch
%%    - mjd
%%    - C1 : Pseudo-range 1 (m)
%%    - C2 : Pseudo-range 2 (m)
%%    - L1 : Phase 1 (cycle number)
%%    - L2 : Phase 2 (cycle number) (for Galileo, L2 corresponds to L5) 
%%  
%% If no informations are founded, no fields are defined in Obs. 
%% Ex : isfield(Obs,'mjd') returns 0 if Obs is empty.
%%  
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Obs_base = struct;
Obs_rover = struct;

% Get corresponding obs in rover RINEX
% 1 - epoch mjd	
mjd_obs = get_mjd_from_epoch(RNX_header_base,epoch);
% 2 -epoch (in rover RINEX) corresponding to mjd_obs
rover_epoch = get_epoch_from_mjd(RNX_header_rover, mjd_obs);
		
if rover_epoch == 0 % no epoch found
	return;
end

% Get observations
Obs_base = get_obs(RNX_header_base,RNX_data_base,constellation,PRN,epoch);		
Obs_rover = get_obs(RNX_header_rover,RNX_data_rover,constellation,PRN,rover_epoch);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
