function run_skyplot(nav,X,Y,Z,mjd_min,mjd_max,cut_off,save)
%% function run_skyplot(nav,X,Y,Z,mjd_min,mjd_max,cut_off,save)
%%
%% Plot a skyplot from orbits
%%
%% Clement Fontaine - 2013-11-19
%%
%% Input : 
%% - nav : nav RINEX file
%% - X, Y, Z : station coordinates
%% - mjd_min, mjd_max : limits of plot
%% - cut_off : cut_off angle (rad)
%% - save : savename (optional)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[NAV_header, NAV_data] = load_rinex_n(nav);

step = (10*60)/86400; % 10min
sat_num = 96; % max number of satellites

Nb_epoch = floor((mjd_max-mjd_min)/step)+1;

mat_skyplot = zeros(sat_num*Nb_epoch,3); % mjd Az El
index = cell(sat_num*Nb_epoch,1); % constPRN

% load data
nb_val = 1;

const = 'GRE';
for const_i = 1:length(const) % constellation
	for PRN = 1:32
	
	tool_print_info(sprintf('Orbit computed : %s%02d',const(const_i),PRN),1);
	
		for mjd = mjd_min:step:mjd_max
		
			% orbits 
			if isfield(NAV_header,'GPSA') % BRDC file
		
				[Eph]=get_ephemeris(NAV_header,NAV_data,const(const_i),PRN,mjd);
				[Xs,Ys,Zs,dte,debug] = orb_sat(Eph,const(const_i),PRN,mjd);
						
			%~ elseif isfield(NAV_header,'Agency') % sp3 file
			
				%~ [X,Y,Z,dte,debug] = orb_sat(NAV_data,const(const_i),PRN,mjd,5);
		
			end
			
			% valid data
			if X~=0
		
				[az,ele,h] = tool_az_ele_h(X,Y,Z,Xs,Ys,Zs);	
				
				if ele>(cut_off*pi/180)
				
					mat_skyplot(nb_val,:) = [mjd,az,ele];
				
					index{nb_val,1} = sprintf('%s%02d',const(const_i),PRN);
						
					nb_val = nb_val + 1;
				
				end

			end
			
		end
	
	end
	
end

%remat
mat_skyplot = mat_skyplot(1:nb_val-1,:);
index = index(1:nb_val-1,:);

% save figure ?
if nargin==9
	if(~strcmp(save,''))
		save = save;
	else
		save = '';
	end
else
	save = '';
end

% skyplot
plot_skyplot(mat_skyplot,index,X,Y,Z,'',save);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
