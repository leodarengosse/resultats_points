function [] = run_plot_obs(rinex_file, const, PRN, save)
%% function [] = run_plot_obs(rinex_file, const, PRN, save)
%%
%% Plot C1, C2, L1, and L2 observations
%%
%% Input : 
%% - rinex_file : name of RINEX obs file
%% - const : constellation ('G' = GPS, 'R' = Glonass, 'E' = Galileo)
%% - PRN : satellite id
%% - save : save name (optional)
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[RNX_header, RNX_data] = load_rinex_o(rinex_file);

nb_ep = RNX_header.nepoch;

mjd = zeros(nb_ep,1);
C1 = zeros(nb_ep,1);
C2 = zeros(nb_ep,1);
L1 = zeros(nb_ep,1);
L2 = zeros(nb_ep,1);

for ep = 1:nb_ep

	obs = get_obs(RNX_header, RNX_data, const, PRN, ep);
	mjd_ep = get_mjd_from_epoch(RNX_header,ep);
	
	if(isfield(obs,'mjd'))
		mjd(ep) = mjd_ep;
		
		if obs.C1~=0
			C1(ep) = obs.C1;
		else
			C1(ep) = NaN;
		end
		
		if obs.C2~=0
			C2(ep) = obs.C2;
		else
			C2(ep) = NaN;
		end
		
		if obs.L1~=0
			L1(ep) = obs.L1;
		else
			L1(ep) = NaN;
		end
		
		if obs.L2~=0
			L2(ep) = obs.L2;
		else
			L2(ep) = NaN;
		end
	else
		
		mjd(ep) = mjd_ep;
		C1(ep) = NaN;
		C2(ep) = NaN;
		L1(ep) = NaN;
		L2(ep) = NaN;
	end
	
end


if(nargin == 4)
	saveC1 = strcat(save,'_C1');
	saveC2 = strcat(save,'_C2');
	saveL1 = strcat(save,'_L1');
	saveL2 = strcat(save,'_L2');
else
	saveC1 = '';
	saveC2 = '';
	saveL1 = '';
	saveL2 = '';
end



figure()
plot_graph(mjd,C1,'C1 pseudo-range','pseudo-range (m)','C1',saveC1);
figure()
plot_graph(mjd,C2,'C2 pseudo-range','pseudo-range (m)','C2',saveC2);
figure()
plot_graph(mjd,L1,'L1 cycle num','cycle number','L1',saveL1);
figure()
plot_graph(mjd,L2,'L2 cycle num','cycle number','L2',saveL2);



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
