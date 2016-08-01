function [] = plot_corr_matrix(Corr,save)
%% function [] = plot_corr_matrix(Corr,save)
%%
%% Print correlation matrix
%%
%% Clement Fontaine 2013-11-14
%%
%% Input : 
%% - Corr : correlation matrix
%% - save : savename (optional)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tool_graphics_toolkit()

if size(Corr,1)==0
	return;
end
%~ figure()

imagesc(abs(Corr));
colormap('gray');

t = colorbar;
ylabel(t, 'Correlation (0 : no correlation, 1 : max correlation)');

% Center image
axis([0.5 size(Corr,1)+0.5 0.5 size(Corr,2)+0.5]);


grid('on')
title('Correlation matrix')

% save as png file
if nargin==2
	if(~strcmp(save,''))
		tool_print_info(sprintf('Correlation matrix saved as : %s',strcat(save,'.png')),1);
		print(strcat(save,'.png'),'-dpng')
	end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
