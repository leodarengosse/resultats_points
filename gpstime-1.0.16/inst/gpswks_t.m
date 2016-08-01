%% Copyright (C) 2010   Jacques Beilin   <jacques.beilin@ensg.eu>
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

function [tgps]=gpswks_t(gpsweek, wsec)
%% -*- texinfo -*-
%% @deftypefn {Function File} {[@var{tgps}] = } gpswks_t (@var{gpsweek},@var{wsec})
%% Calculate all time elements from GPS week and second of week
%%
%% @itemize @bullet
%% @item Input :
%% @itemize @minus
%% @item @var{gpsweek} : GPS week
%% @item @var{wsec} : second of week (0-604800)
%% @end itemize  
%% @item Output :
%% @itemize @minus
%% @item @var{tgps} : gpstime structure
%% @end itemize
%% @end itemize
%%
%% @example
%% example : [tgps]=gpswks_t(1564, 388800)
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

    GPS0 = 315964800.0; % 1980-01-06T00:00:00 in seconds starting from 1970-01-01T00:00:00
    s1970 = gpsweek * 7 * 86400 + wsec + GPS0;

    [tgps]=s1970_t(s1970);

end
