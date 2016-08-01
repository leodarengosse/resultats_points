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

function [tgps]=add_h(tgps,h)
%% -*- texinfo -*-
%% @deftypefn {Function File} {[@var{tgps}] = } add_h (@var{tgps},@var{h})
%% Add @var{h} hours to a tgps structure
%%
%% @itemize @bullet
%% @item Input :
%% @itemize @minus
%% @item @var{tgps} : gpstime structure
%% @item @var{d} : number of hours (+/-)
%% @end itemize  
%% @item Output :
%% @itemize @minus
%% @item @var{tgps} : gpstime structure
%% @end itemize
%% @end itemize
%%
%% @example
%% example : [tgps]=add_h (tgps,1)
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

%% [tgps]=add_h(tgps,h)
%% Adds "h" hours to a gpstime structure
%%
%% Input :
%% - tgps : gpstime structure
%% - h : number of hours (+/-)
%% Output :
%% - tgps : gpstime structure
%%
%% example :
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
%% }

[tgps]=s1970_t(tgps.s1970 + h * 3600);

end
