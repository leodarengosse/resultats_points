function [X,Y,Z] = tool_geocart_GRS80(lon,lat,h)
%% [X,Y,Z] = tool_geocart_GRS80(lon,lat,h)
%%
%% Geographic to cartesian coordinates transformation
%% Jacques BEILIN - 2009-03-19
%% Clement FONTAINE - 2013-10-23
%%
%% Input :
%%	- lon : longitude (decimal degrees)
%%  - lat : latitude (decimal degrees)
%%  - h : height (m)
%%
%% Output : 
%%	- X, Y, Z : cartesian coordinates (m)
%%
%% Function call : 
%% [X,Y,Z] = tool_geocart_GRS80(2.351412386,48.502786872,160.519)
%%
%% Modified :
%% Léo DARENGOSSE - 2016-08-01
%% - Add conversion in radians automatically 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% IAG GRS80 constants
a = 6378137.0;
e2 = 0.006694380022;

%conversion in radians
lat = lat*pi/180; 
lon = lon*pi/180;

% angles in rad
N=a/sqrt(1-e2*(sin(lat))^2);
X=(N+h)*(cos(lon))*(cos(lat));
Y=(N+h)*(sin(lon))*(cos(lat));
Z=(N*(1-e2)+h)*(sin(lat));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
