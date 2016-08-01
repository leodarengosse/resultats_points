function [meanV, stdV] = plot_hist(V,nbin,fig_title,xlab,save);
%% function [meanV, stdV] = plot_hist(V,nbin,fig_title,xlab,save);
%%
%% Plot distribution with a fitted gaussian function
%%
%% Clement Fontaine 2013-11-14
%% 
%% Input : 
%% - V : data vector
%% - nbin : number of bars
%% - fig_title : title
%% - xlab : x label
%% - save : savename (optional)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tool_graphics_toolkit()
if size(V,1)<1
	tool_print_info('No data to plot',3);
	return
end

% New figure
%~ figure()


% mean and standart deviation calculation
meanV = mean(V);
stdV = std(V);

% plot histogram
hist(V,nbin);
[n,xbin]=hist(V,nbin);

% plot gauss function
min_V = min(V);
max_V = max(V);

step = (max_V-min_V)/1000;
x = min_V-3*stdV:step:max_V+3*stdV;

% gauss with normalization
y = (xbin(2)-xbin(1))*length(V)*exp(-(x-meanV).^2./(2*stdV^2))/(stdV*sqrt(2*pi));

hold on
% plot gaussian function
plot(x,y,'r')

%~ % 3 sigma bars
%~ xsigma1 = [meanV - 3*stdV, meanV - 3*stdV];
%~ xsigma2 = [meanV + 3*stdV, meanV + 3*stdV];
%~ ysigma = [0,max(n)];
%~ 
%~ plot(xsigma1,ysigma,'r')
%~ plot(xsigma2,ysigma,'r')

title(fig_title)
xlabel(xlab)
ylabel('Distribution')

%~ legend('Histogram',sprintf('Fitted Gaussian function, mu = %f, sigma = %f',meanV,stdV))
legend('Histogram','Gauss. fit')

% save as png file
if nargin==5
	if(~strcmp(save,''))
	
		tool_print_info(sprintf('Histogram saved as : %s',strcat(save,'.png')),1);
		print(strcat(save,'.png'))
	end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
