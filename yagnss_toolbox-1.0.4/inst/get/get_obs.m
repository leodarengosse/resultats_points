function [Obs] = get_obs(RNX_header,RNX_data,constellation,PRN,epoch)
%% function [Obs] = get_obs(RNX_header,RNX_data,constellation,PRN,epoch) 
%% Get observation data for one satellite at a given epoch (GPS, Glonass and Galileo supported)
%%
%% Clement Fontaine 2013-10-16
%%
%% Input :
%% - RNX_header : structure containing Rinex header
%% - RNX_data : matrix containing data
%%   RNX_header and RNX_data are set up with function load_rinex_o
%% - constellation : 'G' = GPS, 'R' = GLONASS, 'E' = Galileo
%% - PRN : satellite id
%% - epoch
%%
%% Output : 
%% 	- Obs structure containing informations
%%    - constellation
%%    - PRN
%%    - epoch
%%    - mjd
%%    - C1 : Pseudo-range 1 (m)
%%    - C2 : Pseudo-range 2 (m)
%%    - L1 : Phase 1 (cycle number)
%%    - L2 : Phase 2 (cycle number) (for Galileo, L2 corresponds to L5)
%%
%%    
%%  
%% If no informations are founded, no fields are defined in Obs. 
%% Ex : isfield(Obs,'mjd') returns 0 if Obs is empty.
%% 
%%
%% Example :
%%
%% [obs] = get_obs(RNX_header,RNX_data,'G',1,1) for G01
%% obs =
%% {
%%   constellation = G
%%   PRN =  1
%%   epoch =  1
%%   mjd =  56442
%%   C1 =  21325693.1710000 (m)
%%   L1 =  112067331.818000 cycle number
%%   C2 =  21325696.4980000 (m)
%%   L2 =  87325224.7160000 cycle number
%%
%% }
%% 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Output
Obs = struct;

% Test constellation
if strcmp(constellation,'G')
	data = RNX_data.G;
	index = RNX_header.OBS_INDEX_G;
elseif strcmp(constellation,'R')
	data = RNX_data.R;
	index = RNX_header.OBS_INDEX_R;
elseif strcmp(constellation,'E')
	data = RNX_data.E;
	index = RNX_header.OBS_INDEX_E;
else
	tool_print_info('Constellation not implemented : Obs = struct',4);
	return % constellation not implemented
end

[Nsat,Nobs,Nepoch] = size(data);

% Test PRN and epoch number
if ~(PRN<=Nsat && PRN>0 && epoch<=Nepoch && epoch>0)
	tool_print_info('PRN < 1 or PRN > 32 : Obs = struct',4);
	return
end

data_epoch = squeeze(data(PRN,:,epoch)); % RNX_data for one epoch

if(data_epoch(1)>0) % non-0 mjd
	Obs.constellation = constellation;
	Obs.PRN = PRN;
	Obs.epoch = epoch;
	Obs.mjd = data_epoch(1);

	for i=1:size(index,1)
		field = index{i,1};
		Obs.(field) = data_epoch(i+1);		
	end

end

% select data

PR1 = 0;
PR2 = 0;
Phase1 = 0;
Phase2 = 0;


if strcmp(constellation,'G') % GPS
	% RNX v2.11 -> C1, P2, L1, L2
	% RNX v3.02 -> C1C, C2W, L1C, L2W
	
	if isfield(Obs,'C1')
		PR1 = Obs.C1;
	elseif isfield(Obs,'C1C')
		PR1 = Obs.C1C;
	end
	
	if isfield(Obs,'P2')
		PR2 = Obs.P2;
	elseif isfield(Obs,'C2W')
		PR2 = Obs.C2W;
	end
	
	if isfield(Obs,'L1')
		Phase1 = Obs.L1;
	elseif isfield(Obs,'L1C')
		Phase1 = Obs.L1C;
	end	
	
	if isfield(Obs,'L2')
		Phase2 = Obs.L2;
	elseif isfield(Obs,'L2W')
		Phase2 = Obs.L2W;
	end	

elseif strcmp(constellation,'R') % GLONASS
	% RNX v2.11 -> C1, P2, L1, L2
	% RNX v3.02 -> C1C, C2C, L1C, L2C
	
	if isfield(Obs,'C1')
		PR1 = Obs.C1;
	elseif isfield(Obs,'C1C')
		PR1 = Obs.C1C;
	end
	
	if isfield(Obs,'P2')
		PR2 = Obs.P2;
	elseif isfield(Obs,'C2C')
		PR2 = Obs.C2C;
	end

	if isfield(Obs,'L1')
		Phase1 = Obs.L1;
	elseif isfield(Obs,'L1C')
		Phase1 = Obs.L1C;
	end	
	
	if isfield(Obs,'L2')
		Phase2 = Obs.L2;
	elseif isfield(Obs,'L2C')
		Phase2 = Obs.L2C;
	end	


elseif(strcmp(constellation,'E')) % Galileo
	% RNX v2.11 -> C1, C5, L1, L5
	% RNX v3.02 -> C1C or C1X ,C5Q or C5X, L1C or L1X, L5Q or L5X 
	
			
	if isfield(Obs,'C1')
		PR1 = Obs.C1;
	elseif isfield(Obs,'C1C')
		PR1 = Obs.C1C;
	elseif isfield(Obs,'C1X')
		PR1 = Obs.C1X;
	end
	
	if isfield(Obs,'C5')
		PR2 = Obs.C5;
	elseif isfield(Obs,'C5Q')
		PR2 = Obs.C5Q;
	elseif isfield(Obs,'C5X')
		PR2 = Obs.C5X;
	end	
	
		if isfield(Obs,'L1')
		Phase1 = Obs.L1;
	elseif isfield(Obs,'L1C')
		Phase1 = Obs.L1C;
	elseif isfield(Obs,'L1X')
		Phase1 = Obs.L1X;
	end	
	
	if isfield(Obs,'L5')
		Phase2 = Obs.L5;
	elseif isfield(Obs,'L5Q')
		Phase2 = Obs.L5Q;
	elseif isfield(Obs,'L5X')
		Phase2 = Obs.L5X;
	end		


end

Obs.C1 = PR1;
Obs.C2 = PR2;
Obs.L1 = Phase1;
Obs.L2 = Phase2;

if(data_epoch(1)>0) % non-0 mjd
	Obs.C1 = PR1;
	Obs.C2 = PR2;
	Obs.L1 = Phase1;
	Obs.L2 = Phase2;

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
