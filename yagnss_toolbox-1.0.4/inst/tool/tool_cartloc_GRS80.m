function [E,N,U]=  tool_cartloc_GRS80(X0,Y0,Z0,X,Y,Z)
%% function [E,N,U]=  tool_cartloc_GRS80(X0,Y0,Z0,X,Y,Z)
%%
%% Coordinate transformation from cartesian coordinates to local coordinates
%%
%% Clement Fontaine - 2013-12-04
%%
%% Input :
%% - X0, Y0, Z0 : cartesian coordinates of frame origin (m)
%% - X, Y, Z : cartesian coordinates (m)
%%
%% Output : 
%% - E, N, U : Easting, Northing, Up in m
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Geographic coordinates computation
[lon,lat,h] = tool_cartgeo_GRS80(X0,Y0,Z0);

% frame origin
%~ [X1,Y1,Z1] = tool_geocart_GRS80(lon,lat,0);
 
% Local coordinates transformation
R = [ -sin(lon)             cos(lon)            0.0 
      -sin(lat)*cos(lon)   -sin(lat)*sin(lon)   cos(lat)
       cos(lat)*cos(lon)    cos(lat)*sin(lon)   sin(lat) ];
   
        
XYZl = R * [X-X0;Y-Y0;Z-Z0];

E = XYZl(1);
N = XYZl(2);
U = XYZl(3);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
