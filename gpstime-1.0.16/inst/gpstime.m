%% Copyright (C) 2010-2012   Jacques Beilin   <jacques.beilin@ensg.eu>
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

function gpstime
%% Gpstime : GPS date and time manipulation
%% ----------------------------------------
%%
%% gpstime structure example :
%% ---------------------------
%% tgps =
%% {
%%   s1970 =  1262260800        % seconds starting from 1970-01-01T00:00:00
%%   mjd =  55196.0000000000    % Modified Julian Date
%%   wk =  1564                 % GPS week
%%   wsec =  388800             % second of week   
%%   yyyy =  2009               % year   
%%   yy =  9                    % year (2 digits)
%%   mon =  12                  % month
%%   dd =  31                   % day of month
%%   hh =  0                    % hour
%%   min = 0                    % minute
%%   sec = 0                    % second
%%   doy =  365                 % Day of Year
%%   wd =  4                    % Day of week
%%   jd =  2455196.5            % Julian date
%%   dsec =  43200              % second of day
%%   dy =  2009.9993            % decimal year
%% }
%% 
%% Initialization routines : 
%% -------------------------
%%
%% [tgps]=just_now()
%% Calculates all time elements from current computer time
%%
%% [tgps]=gpswkd_t(gpsweek, wday)
%% Calculates all time elements from GPS week and day of week
%%
%% [tgps]=gpswks_t(gpsweek, wsec)
%% Calculates all time elements from GPS week and second of week
%%
%% [tgps]=jd_t(julian_date) 
%% Calculates all time elements from Julian Date
%%
%% [tgps]=mjd_t(mjd) 
%% Calculates all time elements from Modified Julian Date
%%
%% [tgps]=s1970_t(s1970)    
%% Calculates all time elements from seconds starting from 1970-01-01T00:00:00
%%
%% [tgps]=ymdhms_t(y,m,d,hh,mm,ss)
%%  Calculates all time elements from y,m,d,hh,mm,ss
%%
%% [tgps]=yyyyddds_t(yyyy,ddd,s)
%% Calculates all time elements from year, DOY and as an option seconds of day
%%
%% [tgps]=iso_t(iso_time_string)
%% Calculates all time elements from ISO time string (yyyy-mm-ddThh:mm:ss.sss)
%%
%% Addition routines : 
%% -------------------
%%
%% [tgps]=add_s(tgps,s) 
%% Adds "s" seconds to a gpstime structure
%%
%% [tgps]=add_day(tgps,d) 
%% Adds "d" days to a tgps structure
%%
%% Start of period (minute, hour, day, week) calculation
%% -----------------------------------------------------
%%
%% [tgps]=m00(tgps)
%% Fills a tgps structure with current minute at 00
%%
%% [tgps]=h00(tgps)
%% Fills a tgps structure with current hour at 00:00
%%
%% [tgps]=day00(tgps)
%% Fills a tgps structure with current day at 00:00:00
%% 
%% [tgps]=wk00(tgps)
%% Fills a tgps structure with current week at 00:00:00
%%
%% String outputs
%% --------------
%%
%% [s]=st_iso_epoch(tgps,ndec)
%% Outputs ISO time string (yyyy-mm-ddThh:mm:ss.sss)


a=1;
end

