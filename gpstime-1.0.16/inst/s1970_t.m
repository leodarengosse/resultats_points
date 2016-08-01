%% Copyright (C) 2010  Jacques Beilin   <jacques.beilin@ensg.eu>
%%
%% This program is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 3 of the License, or
%% (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, see <http://www.gnu.org/licenses/>
%%
function [tgps]=s1970_t(s1970)
%% -*- texinfo -*-
%% @deftypefn {Function File} {[@var{tgps}] = } s1970_t (@var{s1970})
%% Calculate all time elements from seconds starting from 1970-01-01T00:00:00
%%
%% @itemize @bullet
%% @item Input :
%% @itemize @minus
%% @var{s1970} : number of seconds  from 1970-01-01T00:00:00
%% @end itemize  
%% @item Output :
%% @itemize @minus
%% @item @var{tgps} : gpstime structure
%% @end itemize
%% @end itemize
%% 
%% @example
%% example : [tgps]=s1970_t(1262260800)
%% tgps =
%%    s1970 = 1262260800.00 % seconds starting from 1970-01-01T00:00:00
%%    mjd = 55196.50        % Modified Julian Date
%%    wk = 1564.00          % GPS week
%%    wsec = 388800.00      % second of week 
%%    yyyy = 2009.00        % year (4 digits)
%%    yy = 9.00             % year (2 digits)
%%    mon = 12.00           % month [1-12]
%%    dd = 31.00            % day of month [1-31]
%%    hh = 12.00            % hour
%%    min = 0.00            % minute
%%    sec = 0.00            % second
%%    doy = 365.00          % Day of Year [1-366]
%%    wd = 4.00             % Day of week [0-6]
%%    jd = 2455197.00       % Julian date
%%    jd50 =  21914.50      % Julian date 1950
%%    dsec = 43200.00       % second of day [0-86400]
%%    dy = 2010.00          % decimal year
%%    GMST = 19.34          % Greenwich Mean Sidereal Time
%%    EQEQ = 0.00           % Equation of Equinoxes
%%    GAST = 19.34          % Greenwich Mean Sidereal Time
%% @end example
%%
%% Copyright (C) 2011 IGN/ENSG
%% 
%% Author: Jacques Beilin
%%
%% @end deftypefn



    GPS0 = 315964800.0; % 1980-01-06 en secondes à partir de 1970-01-01
    MJD2000 = 51544.5; % J2000 en mjd
    J2000 = 946728000.0; % J2000 en secondes à partir de 1970-01-01
    J1950 = 2433282.50; % constante de calcul des jours juliens 1950 

    tgps.s1970 = s1970; % secondes à partir de 1970-01-01

    T = gmtime_gps(s1970);

    tgps.mjd = (s1970 - J2000) / 86400 + MJD2000;
    tgps.jd = tgps.mjd + 2400000.5; % jour julien
    tgps.jd50 = tgps.jd - J1950; % JD 1950 (GRGS)

    dsecg = s1970 - GPS0; % secondes depuis GPS0
    njgps = floor ( dsecg / 86400 ) ; % jours entiers depuis GPS0
    tgps.wk = floor ( njgps / 7) ; % semaine GPS
    tgps.wsec = dsecg - tgps.wk *7*86400; % seconde dans la semaine 

    tgps.yyyy = T.year + 1900.0; % annee

    if T.year < 100 
        tgps.yy = tgps.yyyy - 1900;
    else
        tgps.yy = tgps.yyyy - 2000;
    end 
    tgps.mon = T.mon + 1; % mois
    tgps.dd = T.mday; % jour dans le mois

    tgps.hh = T.hour; % heure
    tgps.min = T.min; % minute
    tgps.sec = T.sec + 1e-6 * T.usec; % seconde

    tgps.doy = T.yday +1.0;  % jour de l'année 
    tgps.wd = T.wday; % jour de la semaine (dimanche->0)

    tgps.dsec = (tgps.mjd - floor(tgps.mjd)) * 86400.0; % secondes du jour

    tgps.dy = tgps.yyyy + (tgps.doy - 1  + tgps.dsec / 86400.0) / 365.25; % decimal year

    % Computing GMST
    D = tgps.jd - 2451545.0; % Julian days since J2000
    T = D / 36525.0; % Julian centuries since J2000
    GMST =  18.697374558 + 24.06570982441908 * (tgps.jd - 2451545.0);
    tgps.GMST = mod(GMST,24.0); % Unit : decimal hour

    % Computing Equation of Equinoxes
    d2r = pi / 180.0; 
    Omega = 125.04 - 0.052954 * D;     % Longitude of the ascending node of the Moon
    L = 280.47 + 0.98565 * D;          % Mean Longitude of the Sun
    epsilon = 23.4393 - 0.0000004 * D; % obliquity
    Delta_psi = -0.000319 * sin (Omega*d2r) - 0.000024 * sin (2*L*d2r); % nutation in longitude
    tgps.EQEQ = Delta_psi * cos(epsilon * d2r);

    % Computing GMST
    tgps.GAST = tgps.GMST + tgps.EQEQ;


end
