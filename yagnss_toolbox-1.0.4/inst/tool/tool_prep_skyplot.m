function [mat_skyplot,sat_index] = tool_prep_skyplot(interm_results)
%% function [mat_skyplot,sat_index] = tool_prep_skyplot(interm_results)
%%
%% Set skyplot input from interm_result structure array (see run_spp, run_DGNSS or run_phase for details)
%% 
%% Clement Fontaine
%%
%% Input :
%% - interm_results : structure array containing G and calc, calculated in calc_preprocessing()
%%
%% Output : 
%% - mat_skyplot : matrix contaning [mjd,Az,Elev]
%% - sat_index : cell array containing {'constPRN'} : Format A1I2 : ex : {'G02';'R23'}
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	N_epoch = length(interm_results);
	sat_num = 96;
	
	mat_skyplot = zeros(N_epoch*sat_num,3);
	sat_index = cell(sat_num*N_epoch,2);
	
	N_val = 0;
	for ep = 1:N_epoch
	
		% t, Az and Ele are in result_interm.calc
		nb_loc_val = interm_results(ep).calc.nb_sat;
		
		mat_skyplot(N_val+1:N_val + nb_loc_val,1) = interm_results(ep).calc.t.mjd*ones(nb_loc_val,1);		
		mat_skyplot(N_val+1:N_val + nb_loc_val,2) = interm_results(ep).calc.AzSat(:,1);
		mat_skyplot(N_val+1:N_val + nb_loc_val,3) = interm_results(ep).calc.ElevSat(:,1);
		
		sat_index(N_val+1:N_val + nb_loc_val,1) = interm_results(ep).calc.sat_index(:,1);
		sat_index(N_val+1:N_val + nb_loc_val,2) = interm_results(ep).calc.sat_index(:,2);
		
		N_val = N_val + nb_loc_val;
	
	end
	
	mat_skyplot = mat_skyplot(1:N_val,:);
	sat_index = sat_index(1:N_val,:);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
