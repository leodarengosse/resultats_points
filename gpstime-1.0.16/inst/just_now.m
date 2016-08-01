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

function [tgps]=just_now()
%% [tgps]=just_now() 
%% Calculates all time elements current computer time
%%
%% Input :
%% Output :
%% - tgps : gpstime structure
%%
%% example : [tgps]=just_now() 
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

    s1970 = (now()-719529)*86400; % seconds starting from 1970-01-01T00:00:00

    [tgps]=s1970_t(s1970);
end
