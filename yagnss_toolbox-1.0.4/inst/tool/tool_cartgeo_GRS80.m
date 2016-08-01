function [lon,lat,h] = tool_cartgeo_GRS80(X,Y,Z)
%% function [lon,lat,h] = tool_cartgeo_GRS80(X,Y,Z)
%%
%% Cartesian to geographic coordinates transformation
%%
%% Jacques BEILIN - DPTS - 2012-05-10
%% Clement FONTAINE - 2013-10-23
%%
%% Input : 
%%	- X, Y, Z : cartesian coordinates (m) or vector of coordinates
%%
%% Output : 
%%	- lon : vector of longitude (decimal degrees)
%%	- lat : vector of latitude (decimal degrees)
%%	- h : vector of heights (m)
%%
%% 
%% Function call : 
%% [lamb,phi,h] = tool_cartgeo_GRS80_2([4201575.762;4201575.762],[189856.033;189856.033],[4779066.058;4779066.058])
%%
%% Modified :
%% Léo DARENGOSSE - 2016-08-01
%% - Add conversion in degrees automatically 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% IAG GRS80 constants
a = 6378137.0;
e2 = 0.006694380022;
f=1-sqrt(1-e2);

rxy=sqrt(X.^2+Y.^2);
r=sqrt(X.^2+Y.^2+Z.^2);
mu=atan((Z./rxy).*((1-f)+a.*e2./r));

num=Z.*(1-f)+e2*a.*(sin(mu)).^3;
denum=(1-f)*(rxy-a*e2.*(cos(mu)).^3);
lat=atan(num./denum);

lon = 2 .* atan(Y./(X+rxy));

w=sqrt(1-e2.*sin(lat).^2);
h=rxy.*cos(lat)+Z.*sin(lat)-a.*w;

%%conversion en degrés décimaux
lon = lon*180/pi;
lat = lat*180/pi;
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


