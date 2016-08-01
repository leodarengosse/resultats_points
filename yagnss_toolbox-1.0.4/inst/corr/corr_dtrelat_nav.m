function [dt_relat] = corr_dtrelat_nav(Eph,mjd);
%% [dt_relat] = corr_dtrelat_nav(Eph,mjd);
%% Relativistic correction calculation from navigation message.
%%
%% DPTS - 09/01/2009 - Jacques Beilin
%% Clement Fontaine - 2013-10-23
%%
%% Input 
%% - Eph structure obtained from get_ephemeride_GPS
%%   see get_ephemeride_GPS help for details
%% - mjd : modified Julian date
%%
%% Output
%% - dt_relat : relativistic correction (s)
%%              if no dt_relat computed, 0 is returned
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(isfield(Eph,'alpha0')) % GPS ans Galileo

	% Keplerian elements for satellite k at t
	a = Eph.sqrt_a * Eph.sqrt_a;
	TOC = Eph.TOC; % Time of clock (mjd)
	
	mu = 398600440000000.0;
	n0 = sqrt(mu/a/a/a);
	n = n0 + Eph.delta_n;
	
	% Mean anomaly at t
	tk = (mjd - TOC) * 86400;
	Mk = Eph.M0 + n*tk;
	Ek = Mk;
	while 1==1
	    E1 = Mk + Eph.e * sin(Ek);
	    if (abs(E1-Ek)<1e-12)
	       break   
	    end
	    Ek = E1;
	end
	Ek;
	F = -4.442807633E-10;
	dt_relat = F * Eph.e * sqrt(a) * sin(Ek);
	
elseif (isfield(Eph,'SV_clock_offset')) % Glonass

	%dt_relat = -2 r.v/c^2
	c = 299792458.0;
	[t,Xs,Ys,Zs,VXs,VYs,VZs,dte,debug] = orb_from_RK4(Eph,mjd);
	dt_relat = -2*(Xs*VXs+Ys*VYs+Zs*VZs)/(c*c);
	
else

	tool_print_info('No dt_relat computed',2);
	dte = 0;
	return 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
