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

function [tgps]=m00(tgps)
%% -*- texinfo -*-
%% @deftypefn {Function File} {[@var{tgps}] = } m00 (@var{tgps})
%% Fill a tgps structure with current minute at 00 seconds
%%
%% @itemize @bullet
%% @item Input :
%% @itemize @minus
%% @item @var{tgps} : gpstime structure
%% @end itemize  
%% @item Output :
%% @itemize @minus
%% @item @var{tgps} : gpstime structure
%% @end itemize
%% @end itemize
%%
%% @example
%% example : [tgps]=m00 (tgps)
%% tgps =
%%    s1970 = 1261872000.00 % seconds starting from 1970-01-01T00:00:00
%%    mjd = 55192.00        % Modified Julian Date
%%    wk = 1564.00          % GPS week
%%    wsec = 0.00           % second of week 
%%    yyyy = 2009.00        % year (4 digits)
%%    yy = 9.00             % year (2 digits)
%%    mon = 12.00           % month [1-12]
%%    dd = 27.00            % day of month [1-31]
%%    hh = 0.00             % hour
%%    min = 0.00            % minute
%%    sec = 0.00            % second
%%    doy = 361.00          % Day of Year [1-366]
%%    wd = 4.00             % Day of week [0-6]
%%    jd = 2455192.00       % Julian date
%%    dsec = 43200.00       % second of day [0-86400]
%%    dy = 2009.99          % decimal year
%%    GMST = 6.75           % Greenwich Mean Sidereal Time
%%    EQEQ = 0.00           % Equation of Equinoxes
%%    GAST = 6.75           % Greenwich Mean Sidereal Time
%% @end example
%%
%% Copyright (C) 2011 IGN/ENSG
%% 
%% Author: Jacques Beilin
%% 
%% @end deftypefn

[tgps]=ymdhms_t(tgps.yyyy,tgps.mon,tgps.dd,tgps.hh,tgps.min,0);

end
