function [mjd] = get_mjd_from_epoch(RNX_header,epoch)
%% function [mjd] = get_mjd_from_epoch(RNX_header,epoch)
%% Get mjd from epoch number in RINEX file
%%
%% Clement Fontaine 2013-11-14
%%
%% Input :
%% - RNX_header : Rinex header (see load_rinex_o for details)
%% - epoch : epoch number
%%
%% Output : 
%% - mjd : modified julian day corresponding to epoch
%%  
%%         returns 0 if epoch number not found
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Output
mjd = 0;

if(isfield(RNX_header,'MJD_EPOCH_INDEX'))
	index = RNX_header.MJD_EPOCH_INDEX;
else
	tool_print_info('No field MJD_EPOCH_INDEX in RNX_header : mjd = 0',2);
	return;
end

index_mjd = find(index(:,1)==epoch);

if(length(index_mjd>0))
	mjd = index(index_mjd(1),2);
end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
