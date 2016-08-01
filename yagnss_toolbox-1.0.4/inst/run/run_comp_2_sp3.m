function run_comp_2_sp3(sp31,sp32,sat,mjd_min,mjd_max,step,degree,save)
%% function run_comp_2_sp3(sp31,sp32,sat,mjd_min,mjd_max,step,degree,save)
%%
%% Comparison of sp3 orbits and dte
%%
%% Comparision in NTB (Frenet referential) of second sp3
%%
%% Clement FONTAINE - 2013-11-19
%%
%% Input : 
%% - sp31 : sp3 file 1
%% - sp31 : sp3 file 2
%% - sat : satellite id. Ex 'G23' for GPS satellite of PRN 23
%% - mjd_min, mjd_max : limit for comparison
%% - step : step for orbit computation (in s)
%% - degree for Lagrange interpolation
%% - save : savename (optional)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[sp3_header_1, sp3_data_1] = load_sp3(sp31);
[sp3_header_2, sp3_data_2] = load_sp3(sp32);

if(~isfield(sp3_header_1,'G') || ~isfield(sp3_header_2,'G'))
	return;
end

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

% test degree
if (degree < 1 || degree > 17)
	tool_print_info('1 <= degree <= 17 ',3)
	return
end

step = step/86400; % step in mjd

% Orbit computation
N_epoch = floor((mjd_max-mjd_min)/step)+1;

data = zeros(N_epoch,9); % [mjd,Xsp3_1,Ysp3_1,Zsp3_1,dtesp3_1,Xsp3_2,Ysp3_2,Zsp3_2,dtesp3_2]


% Orbit
nb_val = 1;
for i = mjd_min:step:mjd_max

	[Xsp3_1,Ysp3_1,Zsp3_1,dtesp3_1,debugsp3_1] = orb_sat(sp3_data_1,const,PRN,i,degree);
	[Xsp3_2,Ysp3_2,Zsp3_2,dtesp3_2,debugsp3_2] = orb_sat(sp3_data_2,const,PRN,i,degree);
	
	% save data
	if (Xsp3_1 ==0 || Xsp3_2 == 0)
		data(nb_val,:) = [i,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
	else
		data(nb_val,:) = [i,Xsp3_1,Ysp3_1,Zsp3_1,dtesp3_1,Xsp3_2,Ysp3_2,Zsp3_2,dtesp3_2];
	end
	nb_val = nb_val + 1;
	
end

% redim matrix
data = data(1:nb_val-1,:);

% transformation into inertial referential
% ref time = data(1,1)

data_inert = zeros(size(data));
data_inert(:,1) = data(:,1);

for i = 1:size(data,1)

	theta = (data(i) - data(1))*86400*7.2921151467e-5; % rad, rotation since reference date

	[Xsp3_inert_1,Ysp3_inert_1,Zsp3_inert_1] = tool_rotZ(data(i,2),data(i,3),data(i,4),-theta);
	data_inert(i,2:4) = [Xsp3_inert_1,Ysp3_inert_1,Zsp3_inert_1];
	
	[Xsp3_inert_2,Ysp3_inert_2,Zsp3_inert_2] = tool_rotZ(data(i,6),data(i,7),data(i,8),-theta);
	data_inert(i,6:8) = [Xsp3_inert_2,Ysp3_inert_2,Zsp3_inert_2];
	
	data_inert(i,5) = data(i,5);
	data_inert(i,9) = data(i,9);

end

% local coordinates of sp3 (TBN)
Tr = [data_inert(:,6),data_inert(:,7),data_inert(:,8)];

T = zeros(size(data_inert,1),3);
N = zeros(size(data_inert,1),3);
B = zeros(size(data_inert,1),3);

for i = 1:size(data_inert,1)

	if (i>1 && i<size(data_inert,1))
	
		if (~isnan(data_inert(i,2)) && ~isnan(data_inert(i-1,2))&& ~isnan(data_inert(i+1,2)) && ~isnan(data_inert(i,6)))
		
			%N
	    	N(i,:) = [data_inert(i,6),data_inert(i,7),data_inert(i,8)]./norm([data_inert(i,6),data_inert(i,7),data_inert(i,8)]);
				
			%T
			dxdt = ( data_inert(i+1,6) - data_inert(i-1,6) ) / (2 * ( data_inert(i+1,1) - data_inert(i-1,1) ) * 86400 );
			dydt = ( data_inert(i+1,7) - data_inert(i-1,7) ) / (2 * ( data_inert(i+1,1) - data_inert(i-1,1) ) * 86400 );
			dzdt = ( data_inert(i+1,8) - data_inert(i-1,8) ) / (2 * ( data_inert(i+1,1) - data_inert(i-1,1) ) * 86400 );
			
			dtds = (dxdt^2 + dydt^2 + dzdt^2)^(-0.5);
			
			T(i,1) = dxdt * dtds;
			T(i,2) = dydt * dtds;
			T(i,3) = dzdt * dtds;
	
			%B
			B(i,:) = cross(N(i,:),T(i,:));			
			
		else
		 
			N(i,:) = [NaN,NaN,NaN];
			B(i,:) = [NaN,NaN,NaN];
			T(i,:) = [NaN,NaN,NaN];
			
		end
	else
	
		N(i,:) = [NaN,NaN,NaN];
		B(i,:) = [NaN,NaN,NaN];
		T(i,:) = [NaN,NaN,NaN];
		
	end
	
end

% Transform first SP3 into TBN of second sp3
data_TBN = zeros(size(data_inert,1),4);
for i = 1:size(data_inert,1)

	if(~isnan(N(i,1)))
	
		R = [N(i,:)',T(i,:)',B(i,:)'];
		Tr_i = Tr(i,:)';
		
		Xsp3_1 = data_inert(i,2:4)';
		
		Xsp3_1_sp3_2 = R*(Xsp3_1 - Tr_i);
	
		data_TBN(i,:) = [data_inert(i,1),Xsp3_1_sp3_2'];
		
	else

		data_TBN(i,:) = [data_inert(i,1),NaN,NaN,NaN];
		
	end
end

% save figure ?
if nargin==10
	if(~strcmp(save,''))
		saveX = strcat(save,'X');
		saveY = strcat(save,'Y');
		saveZ = strcat(save,'Z');
		savedte = strcat(save,'dte');
	else
		saveX = '';
		saveY = '';
		saveZ = '';
		savedte = '';
	end
else
	saveX = '';
	saveY = '';
	saveZ = '';
	savedte = '';
end

% plot


% TBN
figure()
subplot(3,1,1)
plot_graph(data_TBN(:,1),data_TBN(:,2),'N_{sp3\_1} - N_{sp3\_2}','d_N (m)','d_N')
subplot(3,1,2)
plot_graph(data_TBN(:,1),data_TBN(:,3),'T_{sp3\_1} - T_{sp3\_2}','d_T (m)','d_T')
subplot(3,1,3)
plot_graph(data_TBN(:,1),data_TBN(:,4),'B_{sp3\_1} - B_{sp3\_2}','d_B (m)','d_B')

if (nargin==8)

	if(~strcmp(save,''))
		tool_print_info(sprintf('Plot saved as : %s',strcat(save,'_pos.png')),1);
		print(strcat(save,'_pos.png'),'-dpng')
	end
	
end



% dte
figure()
plot_graph(data(:,1),[data(:,5)-data(:,9)],'dte_{sp3\_1} - dte_{sp3\_2}','d_{dte} (s)','d_{dte}')
hold on

if (nargin==8)

	if(~strcmp(save,''))
		tool_print_info(sprintf('Plot saved as : %s',strcat(save,'_ddte.png')),1);
		print(strcat(save,'_ddte.png'),'-dpng')
	end
	
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
