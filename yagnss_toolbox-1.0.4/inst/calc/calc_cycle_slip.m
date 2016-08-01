function [t,Pobs,Dobs,PosSat,ElevSat,AzSat,Dtropo,sat_index,amb_index,amb0] = calc_cycle_slip(t,L1,L2,Pobs,Dobs,PosSat,ElevSat,AzSat,Dtropo,sat_index)
%% function [t,Pobs,Dobs,PosSat,ElevSat,Dtropo,sat_index,amb_index,amb0] = calc_cycle_slip(t,L1,L2,Pobs,Dobs,PosSat,ElevSat,Dtropo,sat_index)
%% 
%% Research cycle slips and compute approx values of ambiguities
%%
%% Clement Fontaine 2013-12-17
%%
%% Input : 
%% - t : vector containing mjd corresponding to observations
%% - L1 : matrix containing L1 obs (one column = one receiver)
%% - L2 : matrix containing L2 obs (one column = one receiver)
%% - Pobs : matrix containing L3 obs (one column = one receiver)
%% - Dobs : matrix containing corrected P3 obs (one column = one receiver)
%% - PosSat : matrix containing satellite positions (three columns = one receiver)
%% - ElevSat : matrix containing satellite elevation (one column = one receiver)
%% - AzSat : matrix containing satellite azimuth (one column = one receiver)
%% - Dtropo : matrix containing troposheric delay correction (one column = one receiver)
%% - sat_index : cell array containing satellite id corresponding to observations 
%%   ex : {'G01';'R03';'E15';'G20';'G01';'R03';'E15';'G20'}
%%
%% Output : 
%% Cleaned input data : 
%% - t : vector containing mjd corresponding to observations
%% - L1 : matrix containing L1 obs (one column = one receiver)
%% - L2 : matrix containing L2 obs (one column = one receiver)
%% - Pobs : matrix containing L3 obs (one column = one receiver)
%% - Dobs : matrix containing corrected P3 obs (one column = one receiver)
%% - PosSat : matrix containing satellite positions (three columns = one receiver)
%% - ElevSat : matrix containing satellite elevation (one column = one receiver)
%% - Dtropo : matrix containing troposheric delay correction (one column = one receiver)
%% - sat_index : cell array containing satellite id corresponding to observations
%% Ambiguities : 
%% - amb_index : matrix containing ambiguity index corresponding to observations (one column = one receiver)
%%       Ex : [1 3
%%             2 4
%%			   1 3
%%             2 4
%%			   1 3
%%             2 4];
%% - amb0 : ambiguities approximated values (one line = one ambiguity)
%%       Ex : [12
%%			   1000
%%             30
%%             20] 
%%
%%
%% If less than 15 observations are available to estimate an ambiguity,
%% the function suppress the ambiguity and corresponding obs
%%
%% This function can be used for several receivers (double differences)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nb_obs_per_amb = 20; % 20 obs to estimate an ambiguity
DDL4_lim = 0.0001; % max ampl of DDL4 for cycle slip detection

% cycle split for one or several stations ?
nb_cs = size(L1,2);


amb_index = zeros(size(Pobs,1),nb_cs);
amb0 = zeros(size(Pobs,1),nb_cs);
nb_amb = 0;

t_mjd = unique(t);


for cs = 1:nb_cs

	% L4 = L1 - L2 -- L4 combination
	L4 = L1(:,cs) - L2(:,cs);
	
	obs = [t,L4];
	
	
	const = 'GRE'; % no Glonass -> pour l'instant si
	for const_i = 1:length(const)
	
		for PRN = 1:32
		
			sat = sprintf('%s%02d',const(const_i),PRN);	
			index_pos_sat=find(ismember(sat_index(:,1),sat));	
					
			if(length(index_pos_sat)==0)
				continue
			end
	
			% d2L4/dt2
			obs_sat = obs(index_pos_sat,:); 
			DDL4 = zeros(size(obs_sat,1),1);
			for ep = 2:size(obs_sat,1)-1
				DDL4(ep) = (obs_sat(ep+1,2) - 2*obs_sat(ep,2) + obs_sat(ep-1,2))/((obs_sat(ep+1,1)-obs_sat(ep-1,1))*86400/2)^2;
			end
			
				
			% 1 : obs gap = new ambiguity
			nb_amb = nb_amb + 1;
			amb_index(index_pos_sat(1),cs) = nb_amb; % first obs = new ambiguity
			pos_old_amb = 1;
			
			% first date
			pos_old = find(t_mjd == obs(index_pos_sat(1),1));
			
			% loop on index_pos_sat
			for i = 2:size(index_pos_sat,1)
			
				% gap between obs ?
				pos = find(t_mjd == obs(index_pos_sat(i),1));
				
				if ((pos-pos_old)>1) % gap = new ambiguity
										
					pos_amb = i-1;
					
					%~ tool_print_info(sprintf('Sat %s, amb : %d, min =  %d, max = %d',sat,nb_amb,pos_old_amb,pos_amb),1);
					
					% cycle slip detection
					
					for j = pos_old_amb+1:pos_amb-1
					
					
						% two peacks, opposite signs
						if ( (abs( DDL4(j) ) >= DDL4_lim) && (abs( DDL4(j+1) ) >= DDL4_lim) && ( DDL4(j)/abs(DDL4(j)) + DDL4(j+1)/abs(DDL4(j+1)) == 0 ) )

							cycle_slip_epoch = find(t_mjd == obs_sat(j+1,1));
							tool_print_info(sprintf('Cycle slip detected : Sat %s, epoch : %d',sat,cycle_slip_epoch),1);
							
							nb_amb = nb_amb + 1; % next ambiguity
							amb_index(index_pos_sat(j+1):index_pos_sat(pos_amb),cs) = nb_amb;
							
							
							%~ figure
							%~ plot(DDL4(pos_old_amb:pos_amb))
						
							
						end
						
					
					end
					
					
					
					nb_amb = nb_amb + 1; % next ambiguity
					amb_index(index_pos_sat(i),cs) = nb_amb;
					pos_old_amb = i;
					
				else % no gap
				
					amb_index(index_pos_sat(i),cs) = nb_amb; % same ambiguity
				
				end
				
				pos_old = pos;
				
				
			end	
			
			pos_amb = i;
			
			% if no discontinuities
			for j = pos_old_amb+1:pos_amb-1
					
					
				% two peacks, opposite signs
				if ( (abs( DDL4(j) ) >= DDL4_lim) && (abs( DDL4(j+1) ) >= DDL4_lim) && ( DDL4(j)/abs(DDL4(j)) + DDL4(j+1)/abs(DDL4(j+1)) == 0 ) )

					cycle_slip_epoch = find(t_mjd == obs_sat(j+1,1));
					tool_print_info(sprintf('Cycle slip detected : Sat %s, epoch : %d',sat,cycle_slip_epoch),1);
								
					nb_amb = nb_amb + 1; % next ambiguity
					amb_index(index_pos_sat(j+1):index_pos_sat(pos_amb),cs) = nb_amb;
							
							
					%~ figure
					%~ plot(DDL4(pos_old_amb:pos_amb))
						
							
				end
						
					
			end
			%~ tool_print_info(sprintf('Sat %s, amb : %d, min =  %d, max = %d',sat,nb_amb,pos_old_amb,pos_amb),1);
		
		end
	end		

	
end

% 3 : check if enough obs are available to estimate ambiguity
% remove obs if no enough obs
% amb0 computation


amb_num = 1;
while amb_num <= nb_amb

	[pos_amb,c] = find(amb_index == amb_num);
	
	if length(pos_amb)<nb_obs_per_amb
		L1(pos_amb,:) = [];
		L2(pos_amb,:) = [];
		amb_index(pos_amb,:) = [];
		t(pos_amb) = [];
		Pobs(pos_amb,:) = [];
		Dobs(pos_amb,:) = [];
		PosSat(pos_amb,:) = [];
		ElevSat(pos_amb,:) = [];
		AzSat(pos_amb,:) = [];
		Dtropo(pos_amb,:) = [];
		sat_index(pos_amb,:) = [];
		
		pos_sup = find(amb_index>amb_num);
		amb_index(pos_sup) = amb_index(pos_sup) - 1;
		nb_amb = nb_amb - 1;
	else
	
		amb0(amb_num,:) = mean(Dobs(pos_amb,:) - Pobs(pos_amb,:));
		amb_num = amb_num + 1;
				
	end

end

amb0 = amb0(1:nb_amb,:);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
