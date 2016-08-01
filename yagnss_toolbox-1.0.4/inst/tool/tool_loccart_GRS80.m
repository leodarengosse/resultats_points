function [X,Y,Z]=  tool_loccart_GRS80(X0,Y0,Z0,E,N,U)
%% function [X,Y,Z]=  tool_loccart_GRS80(X0,Y0,Z0,E,N,U)
%%
%% Coordinate transformation from local coordinates to cartesian coordinates
%%
%% Clement Fontaine - 2013-12-19
%%
%% Input :
%% - X0, Y0, Z0 : cartesian coordinates of frame origin (m)
%% - E, N, U : local coordinates (m)
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
   
        
XYZcart = inv(R) * [E; N; U] + [X0; Y0; Z0];

X = XYZcart(1);
Y = XYZcart(2);
Z = XYZcart(3);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
