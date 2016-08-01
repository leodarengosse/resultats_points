function [pivot,interv] = calc_pivot(t,ElevSat,sat_index)
%% function [pivot,interv] = calc_pivot(t,ElevSat,sat_index)
%% Select pivot satellites
%%
%% Clement Fontaine 2013-12-17
%%
%% Input :
%% - t : vector containing mjd corresponding to satellite elevations in ElevSat
%% - ElevSat : vector containing satellite elevation
%% - sat_index : satellite index corresponding to ElevSat data (cell array)
%%
%% Output :
%% - pivot : cell array containing name of pivot satellite
%% - interv : interval of pivot satellite validity [mjd_min mjd_max]
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


pivot = cell(0);
interv = [];

ind = 1;

if size(t,1)==0
    tool_print_info('t is empty : no pivot found',3);
    return
end


% if size(t,1)==0
%     
%     rng = t;
% else
%     
%     rng = t(1):1/24:t(end);
%     
% end

for i = t(1):1/24:t(end)
    
    % pivot satellite -> highest in the sky during 1 hour, present during the whole interval (+- 2 epochs
    
    % obs index between i and i + 1 hour
    obs_id = find((t>=i).*(t<(i+1/24)));
    
    % local obs of current hour
    t_loc = t(obs_id);
    t_unique = unique(t_loc); % mjd without repetitions
    ElevSat_loc = ElevSat(obs_id,:);
    sat_index_loc = sat_index(obs_id,1);
    
    % visible satellites during hour
    satellites = unique(sat_index_loc);
    
    pivot_loc = '';
    Elev_mean_old = 0;
    
    for j = 1:size(satellites,1)
        
        % find ElevSat corresponding to current satellite
        index_ElevSat=find(ismember(sat_index_loc(:,1),satellites(j)));
        Elev_mean = mean(ElevSat_loc(index_ElevSat,1));
        
        % test elevation and visibility
        if (Elev_mean>Elev_mean_old && abs(size(ElevSat_loc(index_ElevSat,1),1)-size(t_unique,1))<=2)
            pivot_loc = satellites(j);
            Elev_mean_old = Elev_mean;
        end
        
    end
    
    tool_print_info(sprintf('Pivot satellite between mjd = %f and mjd = %f : %s',i,i+1/24,pivot_loc{:}),1);
    
    pivot{ind,1} = pivot_loc;
    interv(ind,1:2) = [i,i+1/24];
    
    ind = ind + 1;
    
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
