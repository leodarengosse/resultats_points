function [dte] = corr_dte_nav(Eph,mjd);
%% [dte] = corr_dte_nav(Eph,mjd);
%% Satellite clock error from navigation message.
%%
%% DPTS - 08/05/2012 - Jacques Beilin
%% Clement Fontaine - 2013-10-23
%%
%% Input 
%% - Eph structure obtained from get_ephemeris()
%%   see get_ephemeris help for details
%% - mjd : modified Julian date
%%
%% Output
%% - dte : Satellite clock error (s)
%%
%%         If no dte computed, 0 is returned
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (isfield(Eph,'alpha0')) % GPS and Galileo

	% Difference between mjd and reference orbit date (Time of clock)
	TOC = Eph.mjd; % Time of clock en mjd
	t0 = (mjd - TOC) * 86400;

	% satellite clock error computation
	dte =  Eph.alpha0 + Eph.alpha1 * t0 + Eph.alpha2 * t0^2; 


elseif (isfield(Eph,'SV_clock_offset')) % Glonass

	TOC = Eph.mjd; % Time of clock en mjd
	t0 = (mjd - TOC) * 86400; 
	% satellite clock error computation
	dte = Eph.SV_clock_offset + Eph.SV_relat_freq_offset * t0; % cf doc rinex for signs
else
	tool_print_info('No dte computed',2);
	dte = 0;
	return 
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
