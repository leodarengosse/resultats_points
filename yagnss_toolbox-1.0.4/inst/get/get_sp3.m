function [sp3,nl] = get_sp3(sp3_data,constellation,PRN)
%% function [sp3,nl] = get_sp3(sp3_data,constellation,PRN)
%% Returns sp3 data of a given satellite
%%
%% Clement Fontaine 3013-12-10
%%
%% Input : 
%% - sp3_data : structure created with function load_sp3.m
%% - constellation : 'G' for GPS, 'R' for Glonass and 'E' for Galileo
%% - PRN : satellite id
%%
%% Output : 
%% - sp3 : matrix containing : mjd, X (km), Y (km), Z (km), dte (us)
%% - nl : epoch number
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sp3 = [];
nl = 0;

if ~(strcmp(constellation,'G') || strcmp(constellation,'R') || strcmp(constellation,'E'))
	tool_print_info(sprintf('Constellation %s is not implemented : coordinates and dte set to 0',constellation),2);
	return
end

if(PRN<1 || PRN>32)
	tool_print_info(sprintf('Satellite %s%02d does not exists : coordinates and dte set to 0',constellation,PRN),2);
	return
end

% Get epoch number and sp3 data for selected satellite
if (strcmp(constellation,'G') && isfield(sp3_data,'Gnb'))
	nl = sp3_data.Gnb(PRN);
	sp3 = squeeze(sp3_data.G(PRN,:,:))';
elseif (strcmp(constellation,'R') && isfield(sp3_data,'Rnb'))
	nl = sp3_data.Rnb(PRN);
	sp3 = squeeze(sp3_data.R(PRN,:,:))';
elseif (strcmp(constellation,'E') && isfield(sp3_data,'Enb'))
	nl = sp3_data.Enb(PRN);
	sp3 = squeeze(sp3_data.E(PRN,:,:))';
else
	nl = 0;
end

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
