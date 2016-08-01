function [t,Xs,Ys,Zs,VXs,VYs,VZs,dte,debug]  = orb_from_RK4(Eph,mjd)
%% function [t,Xs,Ys,Zs,VXs,VYs,VZs,dte,debug] = orb_from_RK4(Eph,mjd)
%%
%% Compute position, velocity and dte of a Glonass satellite from its ephemeris at mjd
%%
%% Clement Fontaine 2013-10-29
%%
%% Input :
%% - Eph structure obtained from get_ephemeris
%%   see get_ephemeris help for details
%% - mjd : modified Julian day
%%
%% Output 
%% - t : mjd
%% - Xs,Ys,Zs : cartesian coordinates in WGS84
%% - VXs,VYs,VZs : velocities in WGS84
%% - dte : satellite clock offset
%% - debug : debug structure with all results.
%% 
%% Position, velocity and dte set to 0 if orbit is not computed
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% integration step
tau = 150;
%~ tau = 60;

% 
t=mjd;
Xs = 0;
Ys = 0;
Zs = 0;
VXs = 0;
VYs = 0;
VZs = 0;
dte=0;
debug = cell(0);

% test Eph non vide 
if ~isfield(Eph,'mjd')
	tool_print_info('Eph structure is empty : coordinates and dte set to 0',2);
	return
end

% test fields non 0
if (Eph.X == 0 || Eph.Y == 0 || Eph.Z == 0 || Eph.X_dot == 0 || Eph.Y_dot == 0 || Eph.Z_dot == 0)
	tool_print_info(sprintf('Rinex values = 0 for sat R%02d at %f : coordinates and dte set to 0',Eph.PRN,Eph.mjd),2);
	return
end


% test 'bon' Eph (c a d = GLONASS)
if ~strcmp(Eph.const,'R')
	tool_print_info(sprintf('Bad Constellation : %s, Glonass is required : coordinates and dte set to 0',Eph.const),2);
	return
end

% test Mjd 'pas trop loin' de Eph.mjd
if(abs(mjd-Eph.mjd)>60*30/86400)
	tool_print_info(sprintf('Nearest ephemeris is too far from input mjd : delta = %f s\n',abs(mjd-Eph.mjd)*86400),3);
	return
end

% data

% 
debug.TOE = Eph.mjd;
debug.X = Eph.X;
debug.X_dot = Eph.X_dot;
debug.MS_X_acc = Eph.MS_X_acc;

debug.Y = Eph.Y;
debug.Y_dot = Eph.Y_dot;
debug.MS_Y_acc = Eph.MS_Y_acc;

debug.Z = Eph.Z;
debug.Z_dot = Eph.Z_dot;
debug.MS_Z_acc = Eph.MS_Z_acc;

debug.SV_clock_offset = Eph.SV_clock_offset;
debug.SV_relat_freq_offset = Eph.SV_relat_freq_offset;

% test valeurs ... sur X, Y, Z



% constants
OMEGAe = 7.2921151467e-5; % rad/s
a = 6378136.0;
mu = 398600440000000.0;
C20 = -1.08263e-3;
 
te = mjd_t(debug.TOE);
te_dsec = te.dsec;

% Coordinate transformation to an inertial reference frame

% Sideral time at epoch TOE
debug.theta_Ge = te.GAST*pi/12 + OMEGAe*(te_dsec); % T_rinex = TUTC -> no 3 hours shift

% Position
[debug.Xa,debug.Ya,debug.Za] = tool_rotZ(debug.X,debug.Y,debug.Z,debug.theta_Ge);

% Velocity
[debug.Xa_dot,debug.Ya_dot,debug.Za_dot] = tool_rotZ(debug.X_dot,debug.Y_dot,debug.Z_dot,debug.theta_Ge);
debug.Xa_dot = debug.Xa_dot - OMEGAe*debug.Ya;
debug.Ya_dot = debug.Ya_dot + OMEGAe*debug.Xa;

% Moon and sun acceleration components
[debug.Jms_X,debug.Jms_Y,debug.Jms_Z] = tool_rotZ(debug.MS_X_acc,debug.MS_Y_acc,debug.MS_Z_acc,debug.theta_Ge);

% Runge-Kutta

T = (mjd-te.mjd)*86400;
if(T>0)
	tau = tau;
else
	tau = -tau;
end


N = fix(T/tau) + 2; % 1 for X0 and 1 for rem
rem = (mjd-Eph.mjd)*86400 - fix(T/tau)*tau; % 

%
%    |-----------------------------|
% Eph.mjd                         mjd
%
%    |-----------------------------|
%	 |   |   |   |   |   |   |   | |
%      1   2   3   4   5   6   7  8
%
% 1 to 7 -> tau
% 8 -> rem
%

% matrix initialization

Y0 = zeros(6,1);
Y = zeros(6,N);

K1 = zeros(6,N);
K2 = zeros(6,N);
K3 = zeros(6,N);
K4 = zeros(6,N);

YK = zeros(6,4);

% variable initialization

Y0(1) = debug.Xa;
Y0(2) = debug.Ya;
Y0(3) = debug.Za;

Y0(4) = debug.Xa_dot;
Y0(5) = debug.Ya_dot;
Y0(6) = debug.Za_dot;

Jms_X = debug.Jms_X;
Jms_Y = debug.Jms_Y;
Jms_Z = debug.Jms_Z;

% First column : initial values
Y(:,1) = Y0;
K1(:,1) = Y0;
K2(:,1) = Y0;
K3(:,1) = Y0;
K4(:,1) = Y0;



for i=1:N-1
	
	if(i==N-1)
		tau=rem;
	end

	% K1
	YK(:,1) = Y(:,i);
	
	r = sqrt( YK(1,1)^2 + YK(2,1)^2 + YK(3,1)^2 ) ;
	mu_b = mu / ( r * r ) ;
	x_b = YK(1,1) / r ;
	y_b = YK(2,1) / r ;
	z_b = YK(3,1) / r ;
	rho = a / r ;
	
	K1(:,i+1) = [ YK(4,1)
				  YK(5,1)
				  YK(6,1)
				  - mu_b * x_b + (3/2) * C20 * mu_b * x_b * rho^2 * ( 1 - 5 * z_b ^ 2) + Jms_X
				  - mu_b * y_b + (3/2) * C20 * mu_b * y_b * rho^2 * ( 1 - 5 * z_b ^ 2) + Jms_Y
				  - mu_b * z_b + (3/2) * C20 * mu_b * z_b * rho^2 * ( 3 - 5 * z_b ^ 2) + Jms_Z ];
				  
				  
	% K2
	YK(:,2) = Y(:,i) + 0.5 * tau * K1(:,i+1);
	
	r = sqrt( YK(1,2)^2 + YK(2,2)^2 + YK(3,2)^2 ) ;
	mu_b = mu / ( r * r ) ;
	x_b = YK(1,2) / r ;
	y_b = YK(2,2) / r ;
	z_b = YK(3,2) / r ;
	rho = a / r ;
	
	K2(:,i+1) = [ YK(4,2)
				  YK(5,2)
				  YK(6,2)
				  - mu_b * x_b + (3/2) * C20 * mu_b * x_b * rho^2 * ( 1 - 5 * z_b ^ 2) + Jms_X
				  - mu_b * y_b + (3/2) * C20 * mu_b * y_b * rho^2 * ( 1 - 5 * z_b ^ 2) + Jms_Y
				  - mu_b * z_b + (3/2) * C20 * mu_b * z_b * rho^2 * ( 3 - 5 * z_b ^ 2) + Jms_Z ];
				  	
	% K3
	YK(:,3) = Y(:,i) + 0.5 * tau * K2(:,i+1);
	
	r = sqrt( YK(1,3)^2 + YK(2,3)^2 + YK(3,3)^2 ) ;
	mu_b = mu / ( r * r ) ;
	x_b = YK(1,3) / r ;
	y_b = YK(2,3) / r ;
	z_b = YK(3,3) / r ;
	rho = a / r ;
	
	K3(:,i+1) = [ YK(4,3)
				  YK(5,3)
				  YK(6,3)
				  - mu_b * x_b + (3/2) * C20 * mu_b * x_b * rho^2 * ( 1 - 5 * z_b ^ 2) + Jms_X
				  - mu_b * y_b + (3/2) * C20 * mu_b * y_b * rho^2 * ( 1 - 5 * z_b ^ 2) + Jms_Y
				  - mu_b * z_b + (3/2) * C20 * mu_b * z_b * rho^2 * ( 3 - 5 * z_b ^ 2) + Jms_Z ];
		
	% K4
	YK(:,4) = Y(:,i) + tau * K3(:,i+1);
	
	r = sqrt( YK(1,4)^2 + YK(2,4)^2 + YK(3,4)^2 ) ;
	mu_b = mu / ( r * r ) ;
	x_b = YK(1,4) / r ;
	y_b = YK(2,4) / r ;
	z_b = YK(3,4) / r ;
	rho = a / r ;
	
	K4(:,i+1) = [ YK(4,4)
				  YK(5,4)
				  YK(6,4)
				  - mu_b * x_b + (3/2) * C20 * mu_b * x_b * rho^2 * ( 1 - 5 * z_b ^ 2) + Jms_X
				  - mu_b * y_b + (3/2) * C20 * mu_b * y_b * rho^2 * ( 1 - 5 * z_b ^ 2) + Jms_Y
				  - mu_b * z_b + (3/2) * C20 * mu_b * z_b * rho^2 * ( 3 - 5 * z_b ^ 2) + Jms_Z ];
					  
				  
	% Y(i+1)
	Y(:,i+1) = 	Y(:,i) + (tau / 6) * ( K1(:,i+1) + 2 * K2(:,i+1) + 2 * K3(:,i+1) + K4(:,i+1));		  
	te_dsec = te_dsec + tau;

end

Xs = Y(1,end);
Ys = Y(2,end);
Zs = Y(3,end);
VXs = Y(4,end);
VYs = Y(5,end);
VZs = Y(6,end);


% transformation to PZ-90 (ECEF coordinate system)

debug.theta_Ge_fin = te.GAST*pi/12 + OMEGAe*(te_dsec);

[Xs,Ys,Zs] = tool_rotZ(Xs,Ys,Zs,-debug.theta_Ge_fin);


[VXs,VYs,VZs] = tool_rotZ(VXs,VYs,VZs,-debug.theta_Ge_fin);
VXs = VXs + OMEGAe*Ys;
VYs = VYs - OMEGAe*Xs;

dte = Eph.SV_clock_offset + Eph.SV_relat_freq_offset*(mjd-Eph.mjd)*86400;

% PZ-90.02 to WGS84 (since september 2007 : see navipedia : Reference Frames in GNSS)

T = [ -0.36
       0.08
       0.18 ] ;
       
Res = [Xs;Ys;Zs] + T;

Xs = Res(1);
Ys = Res(2);
Zs = Res(3);

t = floor(te.mjd)+te_dsec/86400;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
