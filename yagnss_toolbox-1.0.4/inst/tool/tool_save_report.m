function [err] = tool_save_report(result,save)
%% function tool_save_report(result,save)
%%
%% Print a result report and information in terminal
%%
%% Clement Fontaine 2013-11-19
%%
%% Input :
%% - result : result structure array (see run_spp run_DGNSS or run_phase for more details)
%% - save : savename
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% open file
[file_id, err ] = fopen(save, 'w');

if(file_id==-1)
    tool_print_info(sprintf('Unable to create : %s',save),3);
    return;
end

if(isempty(result(1).X))
    tool_print_info('No results to display',3);
    return;
end

fprintf(file_id,'---------------------------------------\n');
fprintf(file_id,'                Results                \n');
fprintf(file_id,'---------------------------------------\n');
fprintf(file_id,'\n');
% computation date
fprintf(file_id,'Date : %s\n',datestr(now,'yy_mm_dd_HH:MM:ss'));

fprintf(file_id,'\n');

% results
if(length(result)>1) % normal
    
    tab_pos = zeros(size(result,1),6);
    tab_s02 = zeros(size(result,1),1);
    tab_V = [];
    TO = zeros(size(result,1),2);
    
    fprintf(file_id,'Position and time system offsets \n\n');
    fprintf(file_id,'t(mjd)    X(m)    Y(m)    Z(m)    cdtr(m)    cGGTO(m)    cGPGL(m)\n');
    
    num = 0;
    for i = 1:length(result)
        
        t = result(i).t;
        X = result(i).X;
        Y = result(i).Y;
        Z = result(i).Z;
        cdtr = result(i).cdtr;
        cGGTO = result(i).cGGTO;
        cGPGL = result(i).cGPGL;
        
        if(t==0)
            continue;
        end
        
        num = num + 1;
        tab_pos(num,:) = [X,Y,Z,result(i).E,result(i).N,result(i).U];
        tab_s02(num,1) = result.sigma02;
        tab_V = [tab_V;result(i).V];
        TO(num,:) = [cGGTO,cGPGL];
        
        fprintf(file_id,'%0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f\n',t,X,Y,Z,cdtr,cGGTO,cGPGL);
        
    end
    
    fprintf(file_id,'\n');
    fprintf(file_id,'\n');
    
    fprintf(file_id,'Number of satellites used in compensation \n\n');
    fprintf(file_id,'t(mjd)    nb_GPS    nb_GLO    nb_GAL\n');
    
    for i = 1:length(result)
        t = result(i).t;
        
        if(t==0)
            continue;
        end
        
        nb_GPS = result(i).nb_GPS;
        nb_GLO = result(i).nb_GLO;
        nb_GAL = result(i).nb_GAL;
        
        fprintf(file_id,'%0.3f %d %d %d\n',t,nb_GPS,nb_GLO,nb_GAL);
        
    end
    
    
    fprintf(file_id,'\n');
    fprintf(file_id,'\n');
    
    fprintf(file_id,'DOP quality factors \n\n');
    fprintf(file_id,'t(mjd)    GDOP    PDOP    HDOP    VDOP    TDOP\n');
    
    for i = 1:length(result)
        
        t = result(i).t;
        
        if(t==0)
            continue;
        end
        
        GDOP = result(i).GDOP;
        PDOP = result(i).PDOP;
        HDOP = result(i).HDOP;
        VDOP = result(i).VDOP;
        TDOP = result(i).TDOP;
        
        fprintf(file_id,'%0.3f %0.3f %0.3f %0.3f %0.3f %0.3f\n',t,GDOP,PDOP,HDOP,VDOP,TDOP);
    end
    
    
    % output
    tool_print_info('MEAN POSITION : ',1);
    tool_print_info('',1);
    
    tool_print_info(sprintf('\tX : %0.3f M, std_X : %0.3f M',nanmean(tab_pos(:,1)),nanstd(tab_pos(:,1))),1);
    tool_print_info(sprintf('\tY : %0.3f M, std_Y : %0.3f M',nanmean(tab_pos(:,2)),nanstd(tab_pos(:,2))),1);
    tool_print_info(sprintf('\tZ : %0.3f M, std_Z : %0.3f M',nanmean(tab_pos(:,3)),nanstd(tab_pos(:,3))),1);
    
    tool_print_info('',1);
    
    tool_print_info(sprintf('\tdE : %0.3f M, std_dE : %0.3f M',nanmean(tab_pos(:,4)),nanstd(tab_pos(:,4))),1);
    tool_print_info(sprintf('\tdN : %0.3f M, std_dN : %0.3f M',nanmean(tab_pos(:,5)),nanstd(tab_pos(:,5))),1);
    tool_print_info(sprintf('\tdU : %0.3f M, std_dU : %0.3f M',nanmean(tab_pos(:,6)),nanstd(tab_pos(:,6))),1);
    
    tool_print_info('',1);
    
    tool_print_info('TIME OFFSETS : ',1);
    tool_print_info('',1);
    
    tool_print_info(sprintf('\tcGGTO : %0.3f M, std_cGGTO : %0.3f M',nanmean(TO(:,1)),nanstd(TO(:,1))),1);
    tool_print_info(sprintf('\tcGPGL : %0.3f M, std_cGPGL : %0.3f M',nanmean(TO(:,2)),nanstd(TO(:,2))),1);
    
    tool_print_info('',1);
    
    tool_print_info(sprintf('MEAN SIGMA02 : %0.3f',nanmean(tab_s02)),1);
    
    tool_print_info('',1);
    
    tool_print_info(sprintf('MEAN_RES : %0.3f M',nanmean(tab_V)),1);
    tool_print_info(sprintf('STD_RES : %0.3f M',nanstd(tab_V)),1);
    
    
elseif  (length(result)==1 && size(result.cdtr,1)>1)  % spp pos_glob = 1
    
    fprintf(file_id,'Position and time system offsets \n\n');
    fprintf(file_id,'X(m)    Y(m)    Zm)\n');
    
    X = result.X;
    Y = result.Y;
    Z = result.Z;
    
    fprintf(file_id,'%0.3f %0.3f %0.3f \n',X,Y,Z);
    
    fprintf(file_id,'\n');
    
    fprintf(file_id,'t    cdtr(m)\n');
    
    t = result.t;
    cdtr = result.cdtr;
    
    for i = 1:length(t)
        
        fprintf(file_id,'%0.3f %0.3f \n',t(i),cdtr(i));
        
    end
    
    fprintf(file_id,'\n');
    fprintf(file_id,'cGGTO (m)    cGPGL(m)\n');
    fprintf(file_id,'%0.3f %0.3f \n',result.cGGTO,result.cGPGL);
    
    
    fprintf(file_id,'\n');
    fprintf(file_id,'\n');
    
    fprintf(file_id,'Number of satellites used in compensation \n\n');
    fprintf(file_id,'nb_GPS    nb_GLO    nb_GAL\n');
    
    nb_GPS = result.nb_GPS;
    nb_GLO = result.nb_GLO;
    nb_GAL = result.nb_GAL;
    
    fprintf(file_id,'%d %d %d\n',nb_GPS,nb_GLO,nb_GAL);
    
    
    
    fprintf(file_id,'\n');
    fprintf(file_id,'\n');
    
    fprintf(file_id,'DOP quality factors \n\n');
    fprintf(file_id,'GDOP    PDOP    HDOP    VDOP    TDOP\n');
    
    GDOP = result.GDOP;
    PDOP = result.PDOP;
    HDOP = result.HDOP;
    VDOP = result.VDOP;
    TDOP = result.TDOP;
    
    fprintf(file_id,'%0.3f %0.3f %0.3f %0.3f %0.3f\n',GDOP,PDOP,HDOP,VDOP,TDOP);
    
    % output
    tool_print_info('',1);
    tool_print_info('POSITION : ',1);
    tool_print_info('',1);
    
    tool_print_info(sprintf('\tX : %0.3f M',result.X),1);
    tool_print_info(sprintf('\tY : %0.3f M',result.Y),1);
    tool_print_info(sprintf('\tZ : %0.3f M',result.Z),1);
    
    tool_print_info('',1);
    
    tool_print_info(sprintf('\tdE : %0.3f M',result.E),1);
    tool_print_info(sprintf('\tdN : %0.3f M',result.N),1);
    tool_print_info(sprintf('\tdU : %0.3f M',result.U),1);
    
    tool_print_info('',1);
    
    tool_print_info('TIME OFFSETS : ',1);
    tool_print_info('',1);
    
    tool_print_info(sprintf('\tcGGTO : %0.3f M',result.cGGTO(1)),1);
    tool_print_info(sprintf('\tcGPGL : %0.3f M',result.cGPGL(1)),1);
    
    tool_print_info('',1);
    
    tool_print_info(sprintf('SIGMA02 : %0.3f',result.sigma02),1);
    
    tool_print_info('',1);
    
    tool_print_info(sprintf('MEAN_RES : %0.3f M',nanmean(result.V)),1);
    tool_print_info(sprintf('STD_RES : %0.3f M',nanstd(result.V)),1);
    
    
else
    
    fprintf(file_id,'Rover station \n\n');
    fprintf(file_id,'X_rover(m)    Y_rover(m)    Z_rover(m)\n');
    
    X = result.X;
    Y = result.Y;
    Z = result.Z;
    
    fprintf(file_id,'%0.3f %0.3f %0.3f \n',X,Y,Z);
    
    fprintf(file_id,'\n');
    
    fprintf(file_id,'Number of satellites used in compensation \n\n');
    fprintf(file_id,'nb_GPS    nb_GLO    nb_GAL\n');
    
    nb_GPS = result.nb_GPS;
    nb_GLO = result.nb_GLO;
    nb_GAL = result.nb_GAL;
    
    fprintf(file_id,'%d %d %d\n',nb_GPS,nb_GLO,nb_GAL);
    
    
    
    fprintf(file_id,'\n');
    fprintf(file_id,'\n');
    
    fprintf(file_id,'DOP quality factors \n\n');
    fprintf(file_id,'GDOP    PDOP    HDOP    VDOP    TDOP\n');
    
    GDOP = result.GDOP;
    PDOP = result.PDOP;
    HDOP = result.HDOP;
    VDOP = result.VDOP;
    TDOP = result.TDOP;
    
    fprintf(file_id,'%0.3f %0.3f %0.3f %0.3f %0.3f\n',GDOP,PDOP,HDOP,VDOP,TDOP);
    
    
    % output
    tool_print_info('POSITION : ',1);
    tool_print_info('',1);
    
    tool_print_info(sprintf('\tX : %0.3f M',result.X),1);
    tool_print_info(sprintf('\tY : %0.3f M',result.Y),1);
    tool_print_info(sprintf('\tZ : %0.3f M',result.Z),1);
    
    tool_print_info('',1);
    
    tool_print_info(sprintf('\tdE : %0.3f M, std_E : %0.3f M',result.E, sqrt(result.Qenu(1,1))),1);
    tool_print_info(sprintf('\tdN : %0.3f M, std_E : %0.3f M',result.N, sqrt(result.Qenu(2,2))),1);
    tool_print_info(sprintf('\tdU : %0.3f M, std_E : %0.3f M',result.U, sqrt(result.Qenu(3,3))),1);
    
    tool_print_info('',1);
    
    tool_print_info('TIME OFFSETS : ',1);
    tool_print_info('',1);
    
    tool_print_info(sprintf('\tcGGTO : %0.3f M',result.cGGTO(1)),1);
    tool_print_info(sprintf('\tcGPGL : %0.3f M',result.cGPGL(1)),1);
    
    
    tool_print_info('',1);
    
    tool_print_info(sprintf('SIGMA02 : %0.3f',result.sigma02),1);
    
    tool_print_info('',1);
    
    tool_print_info(sprintf('MEAN_RES : %0.3f M',nanmean(result.V)),1);
    tool_print_info(sprintf('STD_RES : %0.3f M',nanstd(result.V)),1);
    
    
    
    
    
end





fclose(file_id);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%% Nanmean

%% Copyright (C) 2001 Paul Kienzle
%%
%% This program is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 2 of the License, or
%% (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; If not, see <http://www.gnu.org/licenses/>.

%% -*- texinfo -*-
%% @deftypefn {Function File} {@var{v} =} nanmean (@var{X})
%% @deftypefnx{Function File} {@var{v} =} nanmean (@var{X}, @var{dim})
%% Compute the mean value while ignoring NaN values.
%%
%% @code{nanmean} is identical to the @code{mean} function except that NaN values
%% are ignored.  If all values are NaN, the mean is returned as NaN.
%%
%% @seealso{mean, nanmin, nanmax, nansum, nanmedian}
%% @end deftypefn

function v = nanmean (X, varargin)
if nargin < 1
    usage ('v = nanmean(X [, dim])');
else
    n = sum (~isnan(X), varargin{:});
    n(n == 0) = NaN;
    X(isnan(X)) = 0;
    v = sum (X, varargin{:}) ./ n;
end
end

% nanstd


%% Copyright (C) 2001 Paul Kienzle
%%
%% This program is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 2 of the License, or
%% (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; If not, see <http://www.gnu.org/licenses/>.

%% -*- texinfo -*-
%% @deftypefn {Function File} {@var{v} =} nanstd (@var{X})
%% @deftypefnx{Function File} {@var{v} =} nanstd (@var{X}, @var{opt})
%% @deftypefnx{Function File} {@var{v} =} nanstd (@var{X}, @var{opt}, @var{dim})
%% Compute the standard deviation while ignoring NaN values.
%%
%% @code{nanstd} is identical to the @code{std} function except that NaN values are
%% ignored.  If all values are NaN, the standard deviation is returned as NaN.
%% If there is only a single non-NaN value, the deviation is returned as 0.
%%
%% The argument @var{opt} determines the type of normalization to use. Valid values
%% are
%%
%% @table @asis
%% @item 0:
%%   normalizes with @math{N-1}, provides the square root of best unbiased estimator of
%%   the variance [default]
%% @item 1:
%%   normalizes with @math{N}, this provides the square root of the second moment around
%%   the mean
%% @end table
%%
%% The third argument @var{dim} determines the dimension along which the standard
%% deviation is calculated.
%%
%% @seealso{std, nanmin, nanmax, nansum, nanmedian, nanmean}
%% @end deftypefn

function v = nanstd (X, opt, varargin)
if nargin < 1
    usage ('v = nanstd(X [, opt [, dim]])');
else
    if nargin < 3
        dim = min(find(size(X)>1));
        if isempty(dim), dim=1;
        end;
    else
        dim = varargin{1};
    end
    if ((nargin < 2) || isempty(opt))
        opt = 0;
    end
    
    %% determine the number of non-missing points in each data set
    n = sum (~isnan(X), varargin{:});
    
    if n==0
        v = NaN;
        return;
    end
    
    %% replace missing data with zero and compute the mean
    X(isnan(X)) = 0;
    meanX = sum (X, varargin{:}) ./ n;
    
    %% subtract the mean from the data and compute the sum squared
    sz = ones(1,length(size(X)));
    sz(dim) = size(X,dim);
    v = sumsq (X - repmat(meanX,sz), varargin{:});
    
    %% because the missing data was set to zero each missing data
    %% point will contribute (-meanX)^2 to sumsq, so remove these
    v = v - (meanX .^ 2) .* (size(X,dim) - n);
    
    if (opt == 0)
        %% compute the standard deviation from the corrected sumsq using
        %% max(n-1,1) in the denominator so that the std for a single point is 0
        v = sqrt ( v ./ max(n - 1, 1) );
    else if (opt == 1)
            %% compute the standard deviation from the corrected sumsq
            v = sqrt ( v ./ n );
        else
            error ('std: unrecognized normalization type');
        end
        
    end
end
end


function z=sumsq(x,dim)

if(nargin==2)
    z=sum(x.*conj(x),dim);
else
    z=sum(x.*conj(x));
end    

end

