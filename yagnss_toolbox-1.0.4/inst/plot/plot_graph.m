function plot_graph(mjd,data,fig_title,ylab,leg,save,ymin,ymax)
%% function plot_graph(mjd,data,fig_title,ylab,leg,save,ymin,ymax)
%%
%% Plot data
%%
%% Clement Fontaine 2013-11-26
%%
%% 
%% Input :
%% - mjd : mjd vector
%% - data : data matrix
%% - fig_title : title
%% - ylab : y label
%% - leg : legend (comma (,) to set several legends)
%% - save : savename (optional)
%% - ymin, ymax (optional)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tool_graphics_toolkit()
% New figure
%~ figure()

% suppress mjd = 0 if exists
mjd0 = (mjd==0);
mjd(mjd0) = [];
data(mjd0,:) = [];

% test data to plot
if size(data,1)<1
	tool_print_info('No data to plot',3);
	return;
end

%~ % test data variation
%~ if std(data)<0.00001
	%~ data = zeros(size(data));
%~ end

% Offset between local unix time and mjd
gpst = mjd_t(mjd(1));
date = datenum(gpst.yyyy,gpst.mon,gpst.dd,gpst.hh,gpst.min,gpst.sec);
offset = date - mjd(1);

xdate = mjd + offset;

% plot data
hold on
color = ['b','r','g','k','m','c'];
for i = 1:size(data,2)
	plot(xdate,data(:,i),color(floor(i/length(color)) + rem(i,length(color))));
end


% set label
if length(xdate)>=6
	num_ticks = 6;
else
	num_ticks = length(xdate);
end

xdatetick = xdate(1:floor(length(xdate)/num_ticks):end);
labels = datestr((xdatetick),'HH:MM:SS');
set(gca,'XTick',xdatetick)
set(gca,'XTickLabel',labels)

% set ylim if specified
if nargin == 8
	ylim([ymin ymax])
end

title(fig_title)
ylabel(ylab)
grid

leg = cellstr(strsplit(leg,','));
legend(leg)

% save figure
if (nargin == 6 || nargin == 8)
	if(~strcmp(save,''))
		tool_print_info(sprintf('Plot saved as : %s',strcat(save,'.png')),1);
		print(strcat(save,'.png'),'-dpng')
	end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
