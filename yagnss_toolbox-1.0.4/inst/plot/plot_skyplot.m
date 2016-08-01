function [] = plot_skyplot(mat_skyplot,index_sat,Xsta,Ysta,Zsta,sta_name,save)
%% function [] = plot_skyplot(mat_skyplot,index_sat,Xsta,Ysta,Zsta,sta_name,save)
%% Skyplot
%%
%% Clement Fontaine 2013-10-10
%%
%% Input :
%% - mat_skyplot : matrix containing data
%%                 [t,Azimuth,Elevation] : t in mjd, angles in rad
%% - index_sat : cell containing {'constPRN'} conrresponding to
%%               mat_skyplot line. Ex {'G02','R23'}
%% - Xsta, Ysta, Zsta : station coordinates
%% - sta_name : station name
%% - save : savename (optional)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tool_graphics_toolkit()

if length(mat_skyplot)==0
   return;
end
   


figure()
hold on

%~ set(gca, 'defaultTextFontName', 'Arial')

% used constellation
const = 'GRE';

% used colors
color = ['b','r','g'];

% seek Elev = 0
elev_null = find(mat_skyplot(:,2)==0);
mat_skyplot(elev_null,:) = [];
index_sat(elev_null,:) = [];

if size(mat_skyplot,1)<1
	tool_print_info('No data to plot',3);
	return
end

% plot satellites

t_max = max(mat_skyplot(:,1));

for i = 1:length(const)

	for j = 1:32
		
		sat_id = sprintf('%s%02d',const(i),j);	
		index_pos_sat=find(ismember(index_sat(:,1),sat_id));
		
		sat = mat_skyplot(index_pos_sat,:);

		% suppress first if more than three points(if bad station coordinates, first az ele can be far from real position)
		
		if size(sat,1)>3
			sat(1:3,:) = [];
		end

		if(numel(sat)>0)
		
			if size(sat,1)==1
			
				polar(pi/2 - sat(:,2),(pi/2 - sat(:,3))*180/pi,color(i));
			
			else
		
				% cut in order to avoid straight line when satellite reappears
				dt = sat(2:end,1)-sat(1:end-1,1);
				lim = find(dt>1/24);
				lim = [0;lim;size(sat,1)];
				
				for k = 1:length(lim)-1
				
					polar(pi/2 - sat(lim(k)+1:lim(k+1),2),(pi/2 - sat(lim(k)+1:lim(k+1),3))*180/pi,color(i));
				
				end
			
			end
			
			%~ % id at the latest positions, a latest t
			if (sat(end,1) == t_max)
			
				polar(pi/2 - sat(end,2),(pi/2 - sat(end,3))*180/pi,strcat(color(i),'*'));
				[xtxt,ytxt] = pol2cart(pi/2 - sat(end,2),(pi/2 - sat(end,3))*180/pi);
				
				text(xtxt+3,ytxt,sprintf('%s%02d',const(i),j));
			end
			
		end
		
		
	end

end

%%%%%%%%%%%%% CASE OCTAVE %%%%%%%%%%%%%%%%%%%%%%
% plot grid
axis('off')

% Zenithal angle
theta = -pi:2*pi/100:pi;
l = size(theta,2);

for i = 15:15:90

	c = i*ones(1,l);
	polar(theta,c,'k-');
	
	text(0.7*i+2,0.7*i+2,num2str(i))

end

text(-1,100,'0')
text(98,0,'90')
text(-2,-100,'180')
text(-104,0,'370')
%%%%%%%%%%%%% CASE OCTAVE %%%%%%%%%%%%%%%%%%%%%%


axis([-100 100 -100 100])

text(80,90,'G')
plot([85,95],[90,90],'b')

text(80,85,'R')
plot([85,95],[85,85],'r')

text(80,80,'E')
plot([85,95],[80,80],'g')


text(-120,95,sprintf('Station coordinates at t : %0.3f',mat_skyplot(end,1)))
text(-120,90,sprintf('X : %0.3f m',Xsta))
text(-120,85,sprintf('Y : %0.3f m',Ysta))
text(-120,80,sprintf('Z : %0.3f m',Zsta))

title(sprintf('Skyplot (Azimut, Elevation) in Station : %s',sta_name))

% save as png file
if nargin==7
	if(~strcmp(save,''))
		tool_print_info(sprintf('Skyplot saved as : %s',strcat(save,'.png')),1);
		print(strcat(save,'.png'),'-dpng')
	end
end



%%%% TODO -> plot visu sat
%~ figure()
%~ hold on
%~ for i = 1:length(const)
%~ 
	%~ for j = 1:32
		%~ 
		%~ index_pos_sat_temp =(index_const==const(i)).*(index_PRN==j); % 1 if corresponding sat, else 0
		%~ index_pos_sat = find(index_pos_sat_temp==1);
		%~ 
		%~ sat = mat_skyplot(index_pos_sat,:);
%~ 
		%~ if(numel(sat)>0)
		%~ 
			%~ if size(sat,1)==1
			%~ 
				%~ plot(sat(:,1),32*(i-1)+j,'*');
				%~ polar(pi/2 - sat(:,2),(pi/2 - sat(:,3))*180/pi,color(i));
			%~ 
			%~ else
		%~ 
				%~ % cut in order to avoid straight line when satellite reappears
				%~ dt = sat(2:end,1)-sat(1:end-1,1);
				%~ lim = find(dt>1/24);
				%~ lim = [0;lim;size(sat,1)];
				%~ 
				%~ for k = 1:length(lim)-1
				%~ 
					%~ plot(sat(lim(k)+1:lim(k+1),1),(32*(i-1)+j)*ones(length(sat(lim(k)+1:lim(k+1),1)),1),color(i));
					%~ polar(pi/2 - sat(lim(k)+1:lim(k+1),2),(pi/2 - sat(lim(k)+1:lim(k+1),3))*180/pi,color(i));
				%~ 
				%~ end
			%~ 
			%~ end
			%~ 
			%~ % id at the latest positions, a latest t
			%~ if (sat(end,1) == mat_skyplot(end,1))
				%~ polar(pi/2 - sat(end,2),(pi/2 - sat(end,3))*180/pi,strcat(color(i),'*'));
				%~ [xtxt,ytxt] = pol2cart(pi/2 - sat(end,2),(pi/2 - sat(end,3))*180/pi);
				%~ 
				%~ text(xtxt+3,ytxt,sprintf('%s%02d',const(i),j));
			%~ end
			%~ 
		%~ end
		%~ 
		%~ 
	%~ end
%~ 
%~ end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
