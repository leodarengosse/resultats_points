function [] = plot_orbit(NAV_header,NAV_data,sat,mjd_min,mjd_max,save)
%% function [] = plot_orbit(NAV_header,NAV_data,sat,mjd_min,mjd_max,save)
%%
%% Plot ECEF orbits from brdc or sp3
%%
%% Clement Fontaine 2013-11-14
%%
%% Input : 
%% - NAV_header, NAV_data : navigation message (from load_rinex_n)
%%   or sp3_header, sp3_data, precise orbits (from load_sp3)
%% - sat : list of satellites to plot (format : 'R32E12G05')
%% - mjd_min, mjd_max : interval to plot (mjd)
%% - save : savename (optional)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tool_graphics_toolkit()

% earth spinnig velocity
OMEGAe = -7.2921151467e-5; % rad/s


% New figure -> ECEF
ECEF = figure();
hold on
% New figure -> ECI
ECI = figure();
hold on

% read satellite list
if  ~mod (length(sat),3)==0
    tool_print_info('Incorrect satellite list : example :  "G23E12R23"',3);
    return
end

step = (5*60)/86400; % 5 min

% Earth plot
[Xe, Ye, Ze] = sphere (30);
Xe = 6378000*Xe;
Ye = 6378000*Ye;
Ze = 6378000*Ze;
figure(ECEF)
plot3(Xe(:),Ye(:),Ze(:),'k');
figure(ECI)
plot3(Xe(:),Ye(:),Ze(:),'k');

color = ['b','r','g','y','c','m'];
tab_legend = {'Earth'};

% Compute orbit and plot
for i=1:length(sat)/3

	const = sat(3*(i-1)+1);
	PRN = str2num(sat(3*(i-1)+2:3*(i-1)+3));
	
	if ~(strcmp(const,'G') || strcmp(const,'R') || strcmp(const,'E')) 
		continue
	end
	if (PRN>32 || PRN <1)
		continue
	end
	
	
	pos = zeros(floor((mjd_max-mjd_min)/step)+1,7);
	ind = 1;
	
	for j = mjd_min:step:mjd_max
	
		% ECEF
		if isfield(NAV_header,'GPSA') % BRDC file
		
			[Eph]=get_ephemeris(NAV_header,NAV_data,const,PRN,j);
			[X_ECEF,Y_ECEF,Z_ECEF,dte,debug] = orb_sat(Eph,const,PRN,j);
			%~ [X,Y,Z,dte,debug] = pos_sat(Eph,j);
		
		elseif isfield(NAV_header,'Agency') % sp3 file
		
			[Pos_Lagrange]=orb_sp3_Lagrange(NAV_data,j,const,PRN,9);
			X_ECEF = Pos_Lagrange(1);
			Y_ECEF = Pos_Lagrange(2);
			Z_ECEF = Pos_Lagrange(3);
		
		end
		
		% ECI
		[X_ECI,Y_ECI,Z_ECI] = tool_rotZ(X_ECEF,Y_ECEF,Z_ECEF,-OMEGAe*(j-mjd_min)*86400);
		
		pos(ind,:) = [j,X_ECEF,Y_ECEF,Z_ECEF,X_ECI,Y_ECI,Z_ECI];
		ind = ind + 1;

	end
		
	
	% suppress 0
	pos(pos(:,2)==0,:) = [];
	i_col = floor(i/length(color)) + rem(i,length(color));
	
	% plot if satellite exists
	if(length(pos>0))
	
		% ECEF
		figure(ECEF)
		plot3(pos(:,2),pos(:,3),pos(:,4),color(i_col));
		tab_legend(i+1) = {sprintf('%s%02d',const,PRN)};
		
		%ECI
		figure(ECI)
		plot3(pos(:,5),pos(:,6),pos(:,7),color(i_col));
		
	end
	

end

% ECEF
figure(ECEF)

grid

%legend
legend(tab_legend)

%title
title('Satellite positions in ECEF frame')

% axis
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')

axis square

view([pi/4 pi/4 pi/4])

if nargin==6
	if(~strcmp(save,''))
		tool_print_info(sprintf('Orbit plot saved as : %s',strcat(save,'_ECEF.png')),1);
		print(strcat(save,'_ECEF.png'),'-dpng')
	end
end

% ECI
figure(ECI)

grid

%legend
legend(tab_legend)

%title
title('Satellite positions in ECI frame')

% axis
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')

axis square

view([pi/4 pi/4 pi/4])

if nargin==6
	if(~strcmp(save,''))
		tool_print_info(sprintf('Orbit plot saved as : %s',strcat(save,'_ECI.png')),1);
		print(strcat(save,'_ECI.png'),'-dpng')
	end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
