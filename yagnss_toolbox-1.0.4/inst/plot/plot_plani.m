function plot_plani(E,N,leg,save)
%% function plot_plani(E,N,leg,save)
%%
%% Plot E N position
%%
%% Clement Fontaine 2013-11-22
%%
%% Input : 
%% - E : matrix containing E values in column (several columns = several plots)
%% - N : matrix containing N values in column (several columns = several plots)
%% - leg : legend (comma (,) to set several legends)
%% - save : savename (optional)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tool_graphics_toolkit()

hold on

color = ['b','r','g','k','m','c'];

% plot data
for i = 1:size(E,2)
	plot(E(:,i),N(:,i),strcat(color(floor(i/length(color)) + rem(i,length(color))),'.'))
end

% grid
grid

% draw circles
t = linspace(0,2*pi,100)'; 
circsx = cos(t); 
circsy = sin(t); 
plot(circsx,circsy,'r'); 
circsx = 5.*cos(t); 
circsy = 5.*sin(t); 
plot(circsx,circsy,'r'); 
circsx = 10.*cos(t); 
circsy = 10.*sin(t); 
plot(circsx,circsy,'r'); 

% axis limits
axis([-15 15 -15 15]);
axis square

% labels and title
xlabel('E - E0 (m)')
ylabel('N - N0 (m)')
title('dE, dX');

%legend
leg = cellstr(strsplit(leg,','));	
legend(leg)

% save as png file
if nargin==2
	if(~strcmp(save,''))
		tool_print_info(sprintf('Plot2D saved as : %s',strcat(save,'.png')),1);
		print(strcat(save,'.png'),'-dpng')
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
