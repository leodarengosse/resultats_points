function run_comp_sp3(sp3,sat,mjd_min,mjd_max,step,save)
%% function run_comp_sp3(sp3,sat,mjd_min,mjd_max,step,save)
%%
%% Comparison between sp3 orbit (several degrees for Lagrange interpolation)
%%
%% Clement FONTAINE - 2013-12-11
%%
%% Input : 
%% - sp3 : sp3 file
%% - sat : satellite id. Ex 'G23' for GPS satellite of PRN 23
%% - mjd_min, mjd_max : limit for comparison
%% - step : step for orbit computation (in s)
%% - save : savename (optional)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%~ global VERBOSE;
%~ VERBOSE = 1;

[sp3_header, sp3_data] = load_sp3(sp3);

% Test if satellite exists
if( strcmp(sat(1),'G') || strcmp(sat(1),'R') || strcmp(sat(1),'E') )
	if ( str2num(sat(2:3)) > 0 && str2num(sat(2:3)) < 33 )
		const = sat(1);
		PRN = str2num(sat(2:3));
	else
		tool_print_info('Satellite not recognized',3)
		return
	end
else
	tool_print_info('Satellite not recognized',3)
	return
end

step = step/86400; % step in mjd

% Orbit computation
N_epoch = floor((mjd_max-mjd_min)/step)+1;

data = zeros(6, N_epoch,5); % [mjd,Xsp3,Ysp3,Zsp3,dtesp3] fror each degree

% Orbit
degree = [3,5,7,11,13,9];

for d = 1:6
	tool_print_info(sprintf('Computed degree : %d',degree(d)),1)
	nb_val = 1;
	for i = mjd_min:step:mjd_max
					
		% sp3	
		[Xsp3,Ysp3,Zsp3,dtesp3,debugsp3] = orb_sat(sp3_data,const,PRN,i,degree(d));
		
		% save data
		if (Xsp3 == 0)
			data(d,nb_val,:) = [i,NaN,NaN,NaN,NaN];
		else
			data(d,nb_val,:) = [i,Xsp3,Ysp3,Zsp3,dtesp3];
		end
		nb_val = nb_val + 1;
		
	end
end

% save figure ?
if nargin==7
	if(~strcmp(save,''))
		save = save;
	else
		save = '';
	end
else
	save = '';
end

% plot
for i = 1:5

	figure()
	title = sprintf('sp3 comparison : order %d - order 9',degree(i));
	ylab = sprintf('pos_{%d} - pos_{9} (m)',degree(i));
	
	if ~strcmp(save,'');
		savef = sprintf('%s%s%d', save, '_deg', degree(i));
	else
		savef = '';
	end
	
	plot_graph(squeeze(data(1,:,1)),squeeze(data(i,:,2:4))-squeeze(data(6,:,2:4)),title,ylab,'dX, dY, dZ',savef)

end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
