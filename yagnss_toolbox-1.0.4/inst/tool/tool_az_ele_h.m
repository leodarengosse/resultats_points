function [az,ele,h] = tool_az_ele_h(X,Y,Z,Xs,Ys,Zs)
%% fonction [az,ele,h] = tool_az_ele_h(X,Y,Z,Xs,Ys,Zs)
%%
%% Azimuth, elevation and height calculation
%%
%% Jacques BEILIN - ENSG/DPTS - 2014-06-11
%% Clement FONTAINE - 2013-10-23
%%
%% Input :
%%  X,Y,Z : station coordinates
%%  Xs,Ys,Zs : satellite coordinates matrix (1 sat = 1 line)
%% 
%% Output 
%%  az : satellite azimuths vector (rad)
%%  ele : satellite elevation vector (rad)
%%  h : receiver ellipsoidal height (m)
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Geographic coordinates computation
[lon,lat,h] = tool_cartgeo_GRS80(X,Y,Z);

[X1,Y1,Z1] = tool_geocart_GRS80(lon,lat,0);
 
% Local coordinates transformation
R = [ -sin(lon)             cos(lon)            0.0 
      -sin(lat)*cos(lon)   -sin(lat)*sin(lon)   cos(lat)
       cos(lat)*cos(lon)    cos(lat)*sin(lon)   sin(lat) ];
   
az = zeros(length(Xs),1);   
ele = zeros(length(Xs),1);    
        
for i=1:length(Xs)

	XYZl = R * [Xs(i)-X1;Ys(i)-Y1;Zs(i)-Z1];
	D =(XYZl(1)^2+XYZl(2)^2+XYZl(3)^2)^0.5;
	Dh =(XYZl(1)^2+XYZl(2)^2)^0.5;

	% elevation
	E = asin(XYZl(3)/D);
	
	% azimuth
	A = 2 * atan(XYZl(1)/(XYZl(2)+Dh));
	if A < 0 ;	A = A + 2*pi ; end
	
	az(i) = A;
	ele(i) = E;
	
	
	% debug
	%printf('Az = %6.2f , E = %6.2f\n',A *180/pi,E *180/pi);
	
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
