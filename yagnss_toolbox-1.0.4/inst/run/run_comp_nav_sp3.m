function run_comp_nav_sp3(nav,sp3,sat,mjd_min,mjd_max,step,degree,save)
%% function run_comp_nav_sp3(nav,sp3,sat,mjd_min,mjd_max,step,degree,save)
%%
%% Comparison of nav and sp3 orbit and dte
%%
%% Comparision in NTB (Frenet referential)
%%
%% Clement FONTAINE - 2013-11-19
%%
%% Input : 
%% - nav : RINEX nav file
%% - sp3 : sp3 file
%% - sat : satellite id. Ex 'G23' for GPS satellite of PRN 23
%% - mjd_min, mjd_max : limit for comparison
%% - step : step for orbit computation (in s)
%% - degree for Lagrange interpolation
%% - save : savename (optional)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[NAV_header, NAV_data] = load_rinex_n(nav);
[sp3_header, sp3_data] = load_sp3(sp3);

if(~isfield(NAV_header,'TYPE'))
	return;
end

if(~isfield(sp3_header,'G'))
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

data = zeros(N_epoch,9); % [mjd,Xnav,Ynav,Znav,dtenav,Xsp3,Ysp3,Zsp3,dtesp3]

% leap sec (in mjd) if Glonass between nav and sp3
if strcmp(const,'R')
	if(isfield(NAV_header,'LEAP_SECONDS'))
		Leap_sec = NAV_header.LEAP_SECONDS/86400;
	else
		tool_print_info('Field LEAP_SECONDS does not exist',3);
		return;
	end
else
	Leap_sec = 0;
end


%change eph
eph_change = [];
eph_mjd = 0;

% Orbit
nb_val = 1;
for i = mjd_min:step:mjd_max

	% nav
	[Eph]=get_ephemeris(NAV_header,NAV_data,const,PRN,i-Leap_sec);


	[Xnav,Ynav,Znav,dtenav,debugnav] = orb_sat(Eph,const,PRN,i-Leap_sec);
					
	% sp3	
	[Xsp3,Ysp3,Zsp3,dtesp3,debugsp3] = orb_sat(sp3_data,const,PRN,i,degree);
	
	% save data
	if (Xnav ==0 || Xsp3 == 0)
		data(nb_val,:) = [i,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
	else
		data(nb_val,:) = [i,Xnav,Ynav,Znav,dtenav,Xsp3,Ysp3,Zsp3,dtesp3];
	end
	nb_val = nb_val + 1;
	
	if isfield(Eph,'mjd')
		if (~(Eph.mjd==eph_mjd) && (Eph.mjd>mjd_min) && (Eph.mjd<mjd_max))
			eph_change = [eph_change;Eph.mjd];
			eph_mjd = Eph.mjd;
		end
	end
	
end

% redim matrix
data = data(1:nb_val-1,:);

% transformation into inertial referential
% ref time = data(1,1)

data_inert = zeros(size(data));
data_inert(:,1) = data(:,1);

for i = 1:size(data,1)

	theta = (data(i) - data(1))*86400*7.2921151467e-5; % rad, rotation since reference date

	[Xnav_inert,Ynav_inert,Znav_inert] = tool_rotZ(data(i,2),data(i,3),data(i,4),-theta);
	data_inert(i,2:4) = [Xnav_inert,Ynav_inert,Znav_inert];
	[Xsp3_inert,Ysp3_inert,Zsp3_inert] = tool_rotZ(data(i,6),data(i,7),data(i,8),-theta);
	data_inert(i,6:8) = [Xsp3_inert,Ysp3_inert,Zsp3_inert];
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

% Transform brdc into TBN of sp3
data_TBN = zeros(size(data_inert,1),4);
for i = 1:size(data_inert,1)

	if(~isnan(N(i,1)))
	
		R = [N(i,:)',T(i,:)',B(i,:)'];
		Tr_i = Tr(i,:)';
		
		Xbrdc = data_inert(i,2:4)';
		
		Xbrdc_sp3 = R*(Xbrdc - Tr_i);
	
		data_TBN(i,:) = [data_inert(i,1),Xbrdc_sp3'];
		
	else

		data_TBN(i,:) = [data_inert(i,1),NaN,NaN,NaN];
		
	end
end

data_TBN;

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
plot_graph(data_TBN(:,1),data_TBN(:,2),'N_{nav} - N_{sp3}','d_N (m)','d_N')
subplot(3,1,2)
plot_graph(data_TBN(:,1),data_TBN(:,3),'T_{nav} - T_{sp3}','d_T (m)','d_T')
subplot(3,1,3)
plot_graph(data_TBN(:,1),data_TBN(:,4),'B_{nav} - B_{sp3}','d_B (m)','d_B')

if (nargin==8)

	if(~strcmp(save,''))
		tool_print_info(sprintf('Plot saved as : %s',strcat(save,'_pos.png')),1);
		print(strcat(save,'_pos.png'),'-dpng')
	end
	
end



% dte
figure()
plot_graph(data(:,1),[data(:,5)-data(:,9)],'dte_{nav} - dte_{sp3}','d_{dte} (s)','d_{dte}')
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
