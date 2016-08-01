function [] = plot_results_code(result,sta_name,dirsave)
%% function [] = plot_results_code(result,sta_name,dirsave)
%%
%% Plot result from result structure
%% 
%% Clement Fontaine 2013-11-14
%%
%% Input : 
%% - result : array structure containing data to plot 
%% - sta_name : station_name
%% - dirsave : directory to save plots (optional)
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tool_graphics_toolkit()

if length(result)<1
	tool_print_info('No data to plot',3);
	return
end

% Windows or Linux separator
sep = filesep();
if (nargin == 3 && ~strcmp(dirsave,''))
	tool_create_dirs(dirsave);
	dirsave = strcat(dirsave,sep);
else 
	dirsave = '';
end


N_epoch = length(result);
t = zeros(N_epoch,1);
X = zeros(N_epoch,1);
Y = zeros(N_epoch,1);
Z = zeros(N_epoch,1);
cdtr = zeros(N_epoch,1);
cGGTO = zeros(N_epoch,1);
cGPGL = zeros(N_epoch,1);
nb_sat_ac = zeros(N_epoch,3);
Vnorm = [];
DOP = zeros(N_epoch,5);
sigma02 = zeros(N_epoch,1);
E = zeros(N_epoch,1);
N = zeros(N_epoch,1);
U = zeros(N_epoch,1);

sta_pos = result(1).sta_pos;

for i = 1:N_epoch

	% verif si result(i) non vide
	% cas X0 nul... (si pas de lecture du header)
	t(i) = result(i).t;
	X(i) = result(i).X-sta_pos(1);
	Y(i) = result(i).Y-sta_pos(2);
	Z(i) = result(i).Z-sta_pos(3);
	E(i) = result(i).E;
	N(i) = result(i).N;
	U(i) = result(i).U;
	
	cdtr(i) = result(i).cdtr;
	cGGTO(i) = result(i).cGGTO;
	cGPGL(i) = result(i).cGPGL;
	nb_sat(i,:) = [result(i).nb_GPS,result(i).nb_GLO,result(i).nb_GAL];
	Vnorm = [Vnorm;result(i).Vnorm];
	DOP(i,:) = [result(i).GDOP, result(i).PDOP, result(i).HDOP, result(i).VDOP, result(i).TDOP];
	sigma02(i,1) = result(i).sigma02;
		
end

% dirsave 
if (nargin == 3 && ~strcmp(dirsave,''))

	save_XYZ = strcat(dirsave,'XYZ');
	save_ENh = strcat(dirsave,'ENh');
	save_TO = strcat(dirsave,'time_offsets');
	save_quality = strcat(dirsave,'stat_indic');
	save_vis_sat = strcat(dirsave,'vis_sat');
	save_hist_res = strcat(dirsave,'hist_res');
	save_hist_XYZ = strcat(dirsave,'hist_XYZ');
	save_hist_ENh = strcat(dirsave,'hist_ENU');
	save_plani = strcat(dirsave,'pos_plani');



else 
	save_XYZ = '';
	save_ENh = '';
	save_TO = '';
	save_quality = '';
	save_vis_sat = '';
	save_hist_res = '';
	save_hist_XYZ = '';
	save_hist_ENh = '';
	save_plani = '';

end

% position
figure()
subplot(3,1,1)
plot_graph(t,X,sprintf('Station : %s, coordinate variations X - X0 with X0 : %0.3f m',sta_name,sta_pos(1)),'X - X0 (m)','X - X0')
subplot(3,1,2)
plot_graph(t,Y,sprintf('Station : %s, coordinate variations Y - Y0 with Y0 : %0.3f m',sta_name,sta_pos(2)),'Y - Y0 (m)','Y - Y0')
subplot(3,1,3)
plot_graph(t,Z,sprintf('Station : %s, coordinate variations Z - Z0 with Z0 : %0.3f m',sta_name,sta_pos(3)),'Z - Z0 (m)','Z - Z0')

save_plot(save_XYZ);

% E, N, h if defined
if sta_pos(1) ~= 0
	
	figure()
	subplot(3,1,1)
	plot_graph(t,E,sprintf('Station : %s, coordinate variations dE',sta_name),'dE (m)','dE')
	subplot(3,1,2)
	plot_graph(t,N,sprintf('Station : %s, coordinate variations dN',sta_name),'dN (m)','dN')
	subplot(3,1,3)
	plot_graph(t,U,sprintf('Station : %s, coordinate variations dU',sta_name),'dU (m)','dU')
	
	save_plot(save_ENh);
	
	
end

% time offsets
figure()
subplot(3,1,1)
plot_graph(t,cdtr,'c * receiver clock error : cdtr','cdtr (m)','cdtr')
subplot(3,1,2)
plot_graph(t,cGGTO,'c * GPS to Galileo Time Offset : cGGTO','cGGTO (m)','cGGTO')
subplot(3,1,3)
plot_graph(t,cGPGL,'c * GPS to Glonass Time Offset :  cGPGL','cGPGL (m)','cGPGL')

save_plot(save_TO);


% visible satellites before computation and total number of seen satellites
figure()
subplot(2,1,1)
plot_graph(t,nb_sat,'Number of visible satellites per constellation after cut off','Number','GPS,Glonass,Galileo','',min(min(nb_sat))-1,max(max(nb_sat))+1);
subplot(2,1,2)
tot_sat = nb_sat(:,1) + nb_sat(:,2) + nb_sat(:,3);
plot_graph(t,tot_sat,'Number of visible satellites after cut off','Number','Satellite ','',min(min(tot_sat))-1,max(max(tot_sat))+1);

save_plot(save_vis_sat);




% pos histogram
if std(X)>0

	figure()
	subplot(2,2,1)
	plot_hist(X,40,'Histogram : X - X0','Value (m)');
	subplot(2,2,2)
	plot_hist(Y,40,'Histogram : Y - Y0','Value (m)');
	subplot(2,2,3)
	plot_hist(Z,40,'Histogram : Z - Z0','Value (m)');
	subplot(2,2,4)
	axis('off')
	text(0,0.8,sprintf('mean(X-X0) : %0.3f m',mean(X)))
	text(0,0.7,sprintf('std(X-X0) : %0.3f m',std(X)))
	
	text(0,0.5,sprintf('mean(Y-Y0) : %0.3f m',mean(Y)))
	text(0,0.4,sprintf('std(Y-Y0) : %0.3f m',std(Y)))
	
	text(0,0.2,sprintf('mean(Z-Z0) : %0.3f m',mean(Z)))
	text(0,0.1,sprintf('std(Z-Z0) : %0.3f m',std(Z)))


	save_plot(save_hist_XYZ);
	
	if sta_pos(1) ~= 0
	
		figure()
		subplot(2,2,1)	
		plot_hist(E,40,'Histogram : E - E0','Value (m)');
		subplot(2,2,2)	
		plot_hist(N,40,'Histogram : N - N0','Value (m)');
		subplot(2,2,3)	
		plot_hist(U,40,'Histogram : U - U0','Value (m)');
		subplot(2,2,4)
		axis('off')
		text(0,0.8,sprintf('mean(E-E0) : %0.3f m',mean(E)))
		text(0,0.7,sprintf('std(E-E0) : %0.3f m',std(E)))
		
		text(0,0.5,sprintf('mean(N-N0) : %0.3f m',mean(N)))
		text(0,0.4,sprintf('std(N-N0) : %0.3f m',std(N)))
		
		text(0,0.2,sprintf('mean(U-U0) : %0.3f m',mean(U)))
		text(0,0.1,sprintf('std(U-U0) : %0.3f m',std(U)))
	
		save_plot(save_hist_ENh);
		
		
		figure()
		plot_plani(E,N,'Estimated positions');
		save_plot(save_plani)
	
	end

end

% DOP and res histogram
figure()
subplot(2,1,1)
plot_graph(t,DOP,'Dilution of precision','Value','GDOP,PDOP,HDOP,VDOP,TDOP')
subplot(2,1,2)
plot_hist(Vnorm,40,'Histogram : normalized residuals','Value');
%~ plot_graph(t,sigma02(ind_min:end,:),'Square variance factor','','sigma02')

save_plot(save_quality);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function save_plot(save_name)
%% function save_plot(save_name)
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(~strcmp(save_name,''))
	tool_print_info(sprintf('Plot saved as : %s',strcat(save_name,'.png')),1);
	print(strcat(save_name,'.png'),'-dpng')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
