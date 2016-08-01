function [epoch] = get_epoch_from_mjd(RNX_header,mjd,delta)
%% function [epoch] = get_epoch_from_mjd(RNX_header,mjd,delta)
%% Get epoch number from mjd
%%
%% Clement Fontaine 2012-11-13
%%
%% Input :
%% - RNX_header : RNX_header structure (see load_rinex_o for details) 
%% - mjd : modified julian date
%% - delat : max delta if mjd is not exactly found in RNX_header.MJD_EPOCH_INDEX in seconds 
%%           (optional, if not defined = 1e-3)
%%
%% Output : 
%% - epoch : epoch corresponding to mjd
%%           set to 0 is mjd si not found
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Output
epoch = 0;

if(isfield(RNX_header,'MJD_EPOCH_INDEX'))
	index = RNX_header.MJD_EPOCH_INDEX;
else
	tool_print_info('No field MJD_EPOCH_INDEX in RNX_header : epoch = 0',4);
	return;
end

if nargin==2
	delta = 1e-3;
end

% find closest mjd
index_epoch = find(abs(index(:,2)-mjd) == min(abs(index(:,2)-mjd)));

if(abs(index(index_epoch(1),2)-mjd)>delta/86400)
	tool_print_info('Epoch is too far from mjd : epoch = 0',4);
	return
else
	epoch = index(index_epoch(1),1);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
