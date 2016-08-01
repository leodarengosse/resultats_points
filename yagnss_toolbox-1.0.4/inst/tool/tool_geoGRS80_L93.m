function [E,N]= tool_geoGRS80_L93(lambda,phi)
%% function [E,N]=  tool_geoGRS80_L93(lambda,phi)
%% Coordinate transformation from geographic coordinates to Lambert93
%%
%% Clement Fontaine - 2013-11-21
%%
%% Input :
%% - lambda, phi : geographic coordinates in rad (GRS80)
%%
%% Output : 
%% - E, N : Easting, Northing in m
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% IAG GRS 80 (RGF93)
	
	%a = 6378137.0;
	f = 1/298.257222101;

	% first excentricity
	ex = sqrt(2*f-f^2);
		
	% Lambert93 parameters
	lbd0 = 3*pi/180;
	C=11754255.426;
	n=0.7256077650;
	

	
	Xs=700000;
	Ys = 12655612.04987601;

	% isometric latitude
		
	L = log(tan(pi/4+phi/2)*((1-ex*sin(phi))/(1+ex*sin(phi)))^(ex/2));

	R=C*exp(-n*L);
	gam=n*(lambda-lbd0);
	
	% output
	E=Xs+R*sin(gam);
	N=Ys-R*cos(gam);
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

