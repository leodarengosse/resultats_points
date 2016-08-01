function [iono_free_combination] = corr_iono_free(C1,C2,Eph);
%% function [iono_free_combination] = corr_iono_free(C1,C2,Eph);
%% Ionosphere-free combination
%% For GPS, Galileo and Glonass
%%
%% Clement FONTAINE - 2013-10-21
%%
%% Input 
%% 	- C1, C2 : PR 
%%  - Eph : ephemeris
%%
%% Output 
%%	- iono_free : iono_free combination
%%
%%   Returns 0 if const = 'R' and Eph doesn't contain field 'freq_num'
%%   (case of sp3 file)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(C1==0 || C2==0)
	% one obs is not available
	iono_free_combination = 0;
	return;
end


c = 299792458.0;

if strcmp(Eph.const,'G') % GPS L1 and L2
	f1 = 1575.42;
	f2 = 1227.60;
elseif strcmp(Eph.const,'E') % Galileo E1 and E5a
	f1 = 1575.42;
	f2 = 1176.45;
elseif strcmp(Eph.const,'R') % Glonass

	if(isfield(Eph,'freq_num'))
		f1 = 1602.00 + Eph.freq_num*0.5625;
		f2 = 1246.00 + Eph.freq_num*0.4375;
	else
		tool_print_info('Glonass constellation with sp3 file : freq_num is missing, iono_free_combination = 0',2);
		iono_free_combination = 0;
		return
	end
else
	tool_print_info('Constellation not implemented : iono_free_combination = 0',2);
	iono_free_combination = 0;
	return;
end

% gamma 
g = (f1/f2)^2;

iono_free_combination = (C2-g*C1)/(1-g);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
