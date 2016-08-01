function plot_antex(ATX_header, ATX_data, ATX_name, ATX_freq, savename)
%% function plot_antex(ATX_header, ATX_data, ATX_name, ATX_freq, savename)
%%
%% drawing antenna PCV for one antenna
%% 
%% Jacques Beilin - ENSG/DPTS - 2011-12-31
%%
%% Input :
%% - ATX_header : structure containing ANTEX header
%% - ATX_data : cell tab conaining ANTEX data
%% - ATX_name : antenna name (IGS name)
%% - ATX_freq : frequency among ['G01','G02','R01','R02']
%% - savename : save name (optional)
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tool_graphics_toolkit()

% looking for antenna and freq
[ATX] = get_antex(ATX_header, ATX_data, ATX_name, ATX_freq);

if strcmp(ATX.TYPE,'')
	tool_print_info('Structure ATX is empty',3);
	return;
end


d2r = pi / 180;
dzen = (ATX.DAZI(end) - ATX.DAZI(1))/(length(ATX.DAZI)-1);

if dzen>0
	nlignes_azi = size(ATX.VAZI,2);
	m = size(ATX.NOAZI,1); % element number in a NOAZI line

	% conversion in 3 vectors R,THETA,Z
	ATX.R = zeros(nlignes_azi*m,1);
	ATX.THETA = zeros(nlignes_azi*m,1);
	ATX.CORR = zeros(nlignes_azi*m,1);


	for i=1:nlignes_azi
		azi = ATX.VAZI(i);
		for j=1:m
			id = (i-1) * m + j; 
			ATX.R(id) = (j-1) * dzen;
			ATX.THETA(id) = azi * d2r;
			ATX.CORR(id) = ATX.ETAL_AZI(i,j);
		end
	end
			
	% conversion in 3 vectors X,Y,Z
	ATX.X = zeros(nlignes_azi*m,1);
	ATX.Y = zeros(nlignes_azi*m,1);
	ATX.Z = zeros(nlignes_azi*m,1);
	
	for i=1:nlignes_azi
		azi = ATX.VAZI(i);
		for j=1:m
			id = (i-1) * m + j; 
			ATX.X(id) = d2r * ATX.R(id) * cos(ATX.THETA(id));
			ATX.Y(id) = d2r * ATX.R(id) * sin(ATX.THETA(id));
			ATX.Z(id) = ATX.ETAL_AZI(i,j);
		end
	end
			
	step = 4e-2; % grid step
			
	% interpolated grid computation
	xmin=floor(min(ATX.X)/step)*step;
	xmax=ceil(max(ATX.X)/step)*step;
	ymin=floor(min(ATX.Y)/step)*step;
	ymax=ceil(max(ATX.Y)/step)*step;
			
	c = 0;
	nxy = floor((xmax-xmin)*(ymax-ymin)/step/step);
	
	ATX.Xi = zeros(nxy,1);
	ATX.Yi = zeros(nxy,1);
	for x=xmin:step:xmax
		for y=ymin:step:ymax
			c=c+1; 
			ATX.Xi(c) = x;
			ATX.Yi(c) = y;
		end
	end	
	
	ATX.Zi = griddata (ATX.X, ATX.Y, ATX.Z, ATX.Xi, ATX.Yi, 'linear');	
			
	% grid reformat : 2 vectors Xi, Yi and one matrix Zi
	nx = (xmax-xmin)/step +1;
	ny = (ymax-ymin)/step +1;					
	ATX.Xi = zeros(nx,1);
	ATX.Yi = zeros(ny,1);

	c = 0;
	for x=xmin:step:xmax
		c=c+1; 
		ATX.Xi(c) = x;
	end	
		
	c = 0;
	for y=ymin:step:ymax
		c = c + 1;
		ATX.Yi(c) = y;
	end
			
	c=0;
	for x=1:nx
		for y=1:ny
			c=c+1; 
			ATX.ZZi(x,y) = ATX.Zi(c);
		end
	end	
			
	ATX.Zi =ATX.ZZi;
	
	
	figure()
	surf(ATX.Xi,ATX.Yi,ATX.Zi)
	
	% save fig

	if nargin == 5
		if(~strcmp(savename,''))
			tool_print_info(sprintf('Plot saved as : %s',strcat(savename,'_3d.png')),1);
			print(strcat(savename,'_3d.png'),'-dpng')
		end
	end
	
	figure()
	contour(ATX.Xi,ATX.Yi,ATX.Zi)
	
	hold on
	
	% plot grid
	axis('off')

	text(0,2.2,'0')
	text(2.2,0,'90')
	text(-0.1,-2.2,'180')
	text(-2.3,0,'370')
		
		
	% Zenithal angle
	theta = -pi:2*pi/100:pi;
	l = size(theta,2);
	c = 2*ones(1,l);
	polar(theta,c,'k-');
	
	axis([-2.5 2.5 -2.5 2.5])

	if nargin == 5
		if(~strcmp(savename,''))
			tool_print_info(sprintf('Plot saved as : %s',strcat(savename,'_2d.png')),1);
			print(strcat(savename,'_2d.png'),'-dpng')
		end
	end


	end
end	
			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
