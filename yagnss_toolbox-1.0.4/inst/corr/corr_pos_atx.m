function [X, Y, Z] = corr_pos_atx(X, Y, Z, dE, dN, dU, RNX_2_pos)
%% function [X, Y, Z] = corr_pos_atx(X, Y, Z, dE, dN, dU, RNX_2_pos)
%% Phase center and antenna height correction
%%
%% Clement Fontaine - 2014-01-13
%%
%% Input : 
%% - X, Y, Z : position of phase center (if RNX_2_pos == 2) , or station 
%% point (if RNX_2_pos == 1), (m)
%% - dE, dU, dN : corrections (m)
%% - RNX_2_pos : 1 = station point to phase center, 2 = phase center to 
%% station point
%%
%% Output : 
%% - X, Y, Z : position of phase center (if RNX_2_pos == 1) , or station 
%% point (if RNX_2_pos == 2), (m)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global ATX_header;
global ATX_data;

%~ % en var globale ?
antex_file = 'igs08.atx';

if (X~=0)

	% add or supp ENU corrections
	
	if RNX_2_pos == 1
		[X,Y,Z] =  tool_loccart_GRS80(X,Y,Z,dE,dN,dU);
	else
		[X,Y,Z] =  tool_loccart_GRS80(X,Y,Z,-dE,-dN,-dU);
	end

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
