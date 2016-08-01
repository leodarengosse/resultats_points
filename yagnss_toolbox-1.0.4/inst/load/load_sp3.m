function [sp3_header,sp3_data]=load_sp3(filenames)
%% [sp3_Header,sp3_Data]=load_sp3(filename) : sp3 ephemeris loading
%% Interpolation of these data with orb_sp3_Lagrange
%%
%% Jacques Beilin - ENSG/DPTS - 2011-12-17
%% Clement Fontaine - 2013-10-21
%%
%% Input :
%% - filenames : sp3 filenames in cell array(one file for each constellation) : GPS, GLO and GAL supported
%%              ex : {'sp3_GPS.sp3'; 'sp3_GLO.sp3' ; 'sp3_GAL.sp3'}
%%
%% Output :
%% - sp3_header : structure containing sp3 header
%%
%%                   |- Version    
%%                   |- Flag    
%%                   |- Date     
%%                   |- Number_of_Epochs
%%                   |- Data_Used     
%%                   |- Coordinate_Sys    
%%            |- G   |- Orbit_Type  
%% sp3_header |- R   |- Agency
%%            |- E   |- wk  
%%                   |- sow  
%%                   |- Epoch_Interval
%%                   |- mjd
%%                   |- Fractional_Day
%%    
%% - sp3_data : structure containing data
%%  
%%        		  |- G : 3D matrix for GPS satellites
%% sp3_data       |- R : 3D matrix for Glonass satellites
%%                |- E : 3D matrix for Galileo satellites
%%
%%  Matrix index:
%%  	- 1er index : satellite number 
%%      - 2eme index : mjd, X (km), Y (km), Z (km), dte (us)
%%      - 3eme index : several positions, generally one by quarter
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Separator for windows or linux
sep = filesep();

if ischar(filenames)
	filenames = cellstr(filenames);
end

% bin filename
bin_filenames_temp = '';

for i = 1: length(filenames)
	% directory to save bin file
	str = strsplit(filenames{i}(1:end-4),sep);
	dir = '';
	
	for j = 1: length(str)-1
		dir = sprintf('%s%s',dir,str{j},sep);
	end
	
	bin_filenames_temp = sprintf('%s%s_',bin_filenames_temp,str{end});
end

bin_filenames = strcat(dir,bin_filenames_temp(1:end-1),'.sp3.bin');

[Fid,msg]=fopen(bin_filenames);

% debug
 %~ Fid = -1;
% end debug

if Fid>-1
%if Fid==NaN

	filenames_temp = '';
	for i = 1:length(filenames)
		filenames_temp = sprintf('%s %s',filenames_temp,filenames{i});
	end
	
    tool_print_info(sprintf ('Loading binary file  %s (Original sp3 files %s)',bin_filenames,filenames_temp),1);
    fclose(Fid); 
    
    % load
    if exist('OCTAVE_VERSION') % octave 
		load('-binary', bin_filenames, 'sp3_header','sp3_data');
	else
		load (bin_filenames,'-mat', 'sp3_header','sp3_data');
	end

else

	% Output
	sp3_header=cell(0);
	sp3_header.G = cell(0);
	sp3_header.R = cell(0);
	sp3_header.E = cell(0);
	sp3_data=cell(0);
	data.G = [];
	data.R = [];
	data.E = [];
	data.Gnb = [];
	data.Rnb = [];
	data.Enb = [];
	
	for n_fich = 1:length(filenames)
	
	
		sp3_header_temp = cell(0);
		sp3_data_temp = cell(0);

		% File openning
		[fsp3,mess] = fopen(filenames{n_fich}, 'rt');
		% check if file exists
		if (fsp3<=-1)
			tool_print_info(sprintf('Unable to open %s',filenames{n_fich}),3)	
			return;
		    continue; 
		end
	
	
		% header loading
		for i=1:22
			sp3_header_temp.(sprintf('line%02d',i)) = fgetl(fsp3);			
		end
		% file format check
		if  ~strcmp (sp3_header_temp.line01(1:2), '#c') 
		
			tool_print_info('Format not valid : %s',3);
			return; 
		end
	
		tool_print_info(sprintf ('Loading SP3 navigation file %s',filenames{n_fich}),1);
	
		% header reading
		sp3_header_temp.Version = sp3_header_temp.line01(1:2);
		sp3_header_temp.Flag = sp3_header_temp.line01(3:3);
		sp3_header_temp.Date = sp3_header_temp.line01(4:31);
		sp3_header_temp.Number_of_Epochs=sscanf(sp3_header_temp.line01(33:39),'%d');
		
		sp3_header_temp.Data_Used = sp3_header_temp.line01(41:45);
		sp3_header_temp.Coordinate_Sys = sp3_header_temp.line01(47:51);
		sp3_header_temp.Orbit_Type = sp3_header_temp.line01(53:55);
		sp3_header_temp.Agency = sp3_header_temp.line01(57:60);
		
		sp3_header_temp.wk = sscanf(sp3_header_temp.line02(4:7),'%d');
		sp3_header_temp.sow = sscanf(sp3_header_temp.line02(9:23),'%f'); 
		sp3_header_temp.Epoch_Interval = sscanf(sp3_header_temp.line02(25:38),'%f'); 
		sp3_header_temp.mjd = sscanf(sp3_header_temp.line02(40:44),'%d');
		sp3_header_temp.Fractional_Day = sscanf(sp3_header_temp.line02(46:60),'%f'); 
		
		sp3_header_temp.Base_for_Pos_Vel = sscanf(sp3_header_temp.line15(4:13),'%f'); 
		sp3_header_temp.Base_for_Clk_Rate = sscanf(sp3_header_temp.line15(15:26),'%f'); 
	
	
		% data loading
		
		G = zeros(32,5,sp3_header_temp.Number_of_Epochs); % GPS
		Gnb = zeros(32,1);
		
		R = zeros(32,5,sp3_header_temp.Number_of_Epochs); % Glonass
		Rnb = zeros(32,1);
		
		E = zeros(32,5,sp3_header_temp.Number_of_Epochs); % Galileo
		Enb = zeros(32,1);
	
		mjd = 0;
		while ~ feof(fsp3)
			line = fgetl(fsp3);
		    if (strcmp(line(1),'*'))
		    
		        date = str2num(line(2:length(line)));
		        [tgps]=ymdhms_t(date(1),date(2),date(3),date(4),date(5),date(6));
		        mjd = tgps.mjd;
		        
		    elseif (strcmp(line(1:2),'PG'))
			    % GPS   
			    X = sscanf( line(5:18),'%f' );
			    Y = sscanf( line(19:32),'%f' );
			    Z = sscanf( line(33:46),'%f' );
			    dte = sscanf( line(47:60),'%f' );
			    
			    if(~(abs(dte-999999.999999)< 100))
			    
			    	sat = sscanf( line(3:4),'%d' );
					Gnb(sat) = Gnb(sat)+1;
					G(sat,:,Gnb(sat)) =  [ mjd ;  X ; Y ; Z ; dte ];
				end  
		
		    elseif (strcmp(line(1:2),'PR'))
		        % Glonass   
		        X = sscanf( line(5:18),'%f' );
			    Y = sscanf( line(19:32),'%f' );
			    Z = sscanf( line(33:46),'%f' );
			    dte = sscanf( line(47:60),'%f' );
			    
			    if(~(abs(dte-999999.999999)< 100))
					sat = sscanf( line(3:4),'%d' );
			        Rnb(sat) = Rnb(sat)+1;  
			       	R(sat,:,Rnb(sat)) =  [ mjd ;  X ; Y ; Z ; dte ]; 
				end
				
		    elseif (strcmp(line(1:2),'PE'))
		        % Galileo
		        X = sscanf( line(5:18),'%f' );
			    Y = sscanf( line(19:32),'%f' );
			    Z = sscanf( line(33:46),'%f' );
			    dte = sscanf( line(47:60),'%f' );
			   
			    if(~(abs(dte-999999.999999)< 100))
			        sat = sscanf( line(3:4),'%d' );
			        Enb(sat) = Enb(sat)+1; 	
			       	E(sat,:,Enb(sat)) =  [ mjd ;  X ; Y ; Z ; dte ];
				end
				
		    elseif (strcmp(line(1:3),'EOF'))
		    end%if   
		end%while
		
		% file closing
		sp3_header_temp.name=filenames{n_fich};
		status = fclose(fsp3);
		
		if(sum(Gnb>0))
			sp3_data.G = G;
			sp3_data.Gnb = Gnb;
			sp3_header.G = sp3_header_temp;
		end
		if(sum(Rnb>0))
			sp3_data.R = R;
			sp3_data.Rnb = Rnb;
			sp3_header.R = sp3_header_temp;
		end
		if(sum(Enb>0))
			sp3_data.E = E;
			sp3_data.Enb = Enb;
			sp3_header.E = sp3_header_temp;
		end
		

		
	end
	
	% save
	if exist('OCTAVE_VERSION') % octave 
		save('-binary', bin_filenames, 'sp3_header','sp3_data');
	else
		save(bin_filenames,'-mat', 'sp3_header','sp3_data');      
	end


end 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

