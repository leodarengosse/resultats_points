function [RNX_header,RNX_data]=load_rinex_o(filename,epoch_max)
%% function [RNX_header,RNX_data]=load_rinex_o(filename,epoch_max)
%% GPS/Glonass/Galileo observation RINEX loading
%% 
%% Jacques Beilin - ENSG/DPTS - 2011-12-21
%% Clement Fontaine - 2013-10-14
%%
%% RINEX prereqisites :
%% - RINEX v2.__ or v3.__
%% - only 1 header bloc
%% - GPS,GLONASS and Galileo. Other systems will be ignored.
%% - only epoch with code=0 or 1 will be loaded. Other epochs will be ignored.
%%
%% Input 
%%
%% - filename : RINEX file name
%% - (epoch_max) : OPTIONAL, number epoch to be loaded
%% 
%% Output
%%
%% - RNX_header = structure containing RINEX header elements
%%   {
%%  VERSION =  3.02000000000000
%%  TYPE = M
%%  MARKER_NAME = VILL
%%  MARKER_NUMBER = 13406M001
%%  REC_N = 3001316             
%%  REC_TYPE = SEPT POLARX4        
%%  ANT_N = 5166                
%%  ANT_TYPE = SEPCHOKE_MC     NONE
%%  X =  4201792.29500000
%%  Y =  177945.238000000
%%  Z =  4779286.68500000
%%  dH =  0.0937000000000000
%%  dE = 0
%%  dN = 0
%%  INTERVAL =  30
%%  MJD_EPOCH_INDEX = [mjd,epoch]
%%		56442 1
%%      56442.01 2
%%  OBS_INDEX_G,OBS_INDEX_R,OBS_INDEX_E =
%%  	Structures with : {'observation type',index for this constellation} f
%%      for GPS, Galileo and Glonass
%%
%%      Corresponds to field order in RINEX file.
%%  nepoch =  15
%%
%% - RNX_data : structure containing data (RNX_data.G = GPS, RNX_data.R = Glonass, RNX_data.E = Galileo)
%%   (3D matrix 32 satellites x n observables x 3600 epochs)
%%   All satellites which are not observed will be set to 0
%%   Index order : 
%%   - 1st index : satellite PRN  
%%   - 2nd index : mjd + observables
%%     Date is provided for each PRN while commpon for all satellite. It can be used to 
%%     find observed satellites.
%%   - 3rd index : epoch. Epochs are added without any time considerations.
%%
%%	INFO : 
%%	Use the function get_obs() to find data in RNX_dataX

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bin_filename = sprintf('%s.bin',filename);

[Fid,msg]=fopen(bin_filename);


% debug 
 %~ Fid=-1; % do not read .bin even if it exists
% end debug 

if Fid>-1
%if Fid==NaN

    tool_print_info(sprintf ('Loading binary file  %s (Original RINEX file %s)',bin_filename,filename),1);
    fclose(Fid); 
    
    % load
    if exist('OCTAVE_VERSION') % octave 
		load('-binary', bin_filename, 'RNX_header','RNX_data');
	else % matlab
		load (bin_filename,'-mat', 'RNX_header','RNX_data');
	end
        
    % test epoch ?

else

    % initialization of variables
    RNX_header=[];            % RINEX header
    RNX_header.VERSION=0;  
    RNX_data=[];
    RNX_dataG=[];             % GPS data
    RNX_dataR=[];             % Glonass data
    RNX_dataE=[];             % Galileo data
    
    

    % opening of RINEX file 
    [frnx,mess] = fopen(filename, 'rt');

    % check if file exists
    if frnx<=-1
		tool_print_info(sprintf('Unable to open %s',filename),3)	
		return;
    else
    	tool_print_info(sprintf('Loading RINEX file %s',filename),1);

        % header loading
        [RNX_header]=read_header(frnx,RNX_header);

        % data loading        
        if nargin == 1
			epoch_max = -1;
		end
        
        if RNX_header.VERSION>=3.0
			[RNX_dataG,RNX_dataR,RNX_dataE,RNX_header.nepoch]=read_data_3(frnx,RNX_header,epoch_max); % RINEX v3.__
        else
			[RNX_dataG,RNX_dataR,RNX_dataE,RNX_header.nepoch]=read_data_2_11(frnx,RNX_header,epoch_max); % Rinex v2.11
        end
        
        RNX_data.G = RNX_dataG;
        RNX_data.R = RNX_dataR;
        RNX_data.E = RNX_dataE;
        
        % clean data
        [RNX_data] = qc(RNX_data);
        
        % Set mjd_epoch_index
        mjd_epoch_index = zeros(size(RNX_data.G,3),2);
       
        
        for i = 1:size(RNX_data.G,3)
        
			mjd_loc = [max(RNX_data.G(:,1,i)) max(RNX_data.R(:,1,i)) max(RNX_data.E(:,1,i))];
			mjd_epoch_index(i,:) = [i, max(mjd_loc)];
        
        end
        
        RNX_header.MJD_EPOCH_INDEX = mjd_epoch_index;
        
        status = fclose(frnx);
        
        % save binary
		if exist('OCTAVE_VERSION') % octave 
			save('-binary', bin_filename, 'RNX_header','RNX_data');
		else % matlab
			save(bin_filename,'-mat', 'RNX_header','RNX_data');      
		end
            
    end
      
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [RNX]=read_header(frnx,RNX)
%% function [RNX]=read_header(frnx,RNX) : RINEX header loading
%% 
%% Jacques Beilin - ENSG/DPTS - 2011-12-18
%% Clement Fontaine - 2013-10-14
%%
%% Input :
%% - frnx : rinex file handle (already openned)
%% - RNX : Rinex header structure
%%
%% Output :
%% - RNX : Rinex header structure
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RNX.VERSION=0;
RNX.TYPE='';
RNX.MARKER_NAME='UNKN';
RNX.MARKER_NUMBER='';
RNX.REC_N='';
RNX.REC_TYPE='';
RNX.ANT_N='';
RNX.ANT_TYPE='';
RNX.X=0;
RNX.Y=0;
RNX.Z=0;
RNX.dH=0;
RNX.dE=0;
RNX.dN=0;
RNX.INTERVAL=0;
RNX.OBS_INDEX_G=cell(0);
RNX.OBS_INDEX_R=cell(0);
RNX.OBS_INDEX_E=cell(0);
RNX.MJD_EPOCH_INDEX = [];
	
constellation = ''; 
    
i=1;
while (i<100) 
	line=fgetl(frnx);
	header(i,1:length(line)) = line;
	i=i+1;
	%printf('>%s<\n',header(i-1,61:73))
	
	% RINEX version number and type
	if strcmp(header(i-1,61:80),   'RINEX VERSION / TYPE') 
		RNX.VERSION=str2double(header(i-1,1:10));
		RNX.TYPE=header(i-1,41:41);
	end
	
	% Station name
	if strcmp(header(i-1,61:71),   'MARKER NAME') ; RNX.MARKER_NAME=header(i-1,1:4) ; end
	
	% Station number
	if strcmp(header(i-1,61:73),   'MARKER NUMBER') ; RNX.MARKER_NUMBER=header(i-1,1:9) ; end
	
	% Receiver number and type
	if strcmp(header(i-1,61:79),   'REC # / TYPE / VERS') 
		RNX.REC_N=header(i-1,1:20);
		RNX.REC_TYPE=header(i-1,21:40);
	end
	
	% Antenna type
	if strcmp(header(i-1,61:72),   'ANT # / TYPE') 
		RNX.ANT_N=header(i-1,1:20);
	   RNX.ANT_TYPE=header(i-1,21:40);
	end
	
	% Approximated position
	if strcmp(header(i-1,61:79),   'APPROX POSITION XYZ') 
	
		% Official Rinex Format
		%RNX.X=str2double(header(i-1,1:14));
		%RNX.Y=str2double(header(i-1,15:28));
		%RNX.Z=str2double(header(i-1,29:42));
		
		% ... this works for any Format (pb with mgex data)
		[val,count] = sscanf(header(i-1,1:60),'%f %f %f');
		RNX.X=val(1);
		RNX.Y=val(2);
		RNX.Z=val(3);
		
	end
	
	% Antenna : delta
	if strcmp(header(i-1,61:80),   'ANTENNA: DELTA H/E/N') 
	
		% Official Rinex Format
		%RNX.dH=str2double(header(i-1,1:14));
		%RNX.dE=str2double(header(i-1,15:28));
		%RNX.dN=str2double(header(i-1,29:42));
		
		% ... this works for any Format (pb with mgex data)
		[val,count] = sscanf(header(i-1,1:60),'%f %f %f');
		RNX.dH=val(1);
		RNX.dE=val(2);
		RNX.dN=val(3);
		
	end
	
	% Interval
	if strcmp(header(i-1,61:68),   'INTERVAL') 
				
		RNX.INTERVAL = str2double(header(i-1,1:10));
		
	end
	
	% List of observable
	if strcmp(header(i-1,61:79),   '# / TYPES OF OBSERV') % Case of RINEX V2.11
		temp = regexprep(strtrim(header(i-1,11:60)),'\s{1,10}',' ');
		tab_obs = strsplit(temp,' ');
		
		% Filling of RNX_OBS index	
		for j=1:size(tab_obs,2)
		
			% idem for GPS, GLONASS and Galileo
			ind_obs = size(RNX.OBS_INDEX_G,1)+1;
			RNX.OBS_INDEX_G{ind_obs,1} = tab_obs{1,j};
			RNX.OBS_INDEX_G{ind_obs,2} = ind_obs;
			
			ind_obs = size(RNX.OBS_INDEX_R,1)+1;
			RNX.OBS_INDEX_R{ind_obs,1} = tab_obs{1,j};
			RNX.OBS_INDEX_R{ind_obs,2} = ind_obs;
			
			ind_obs = size(RNX.OBS_INDEX_E,1)+1;
			RNX.OBS_INDEX_E{ind_obs,1} = tab_obs{1,j};
			RNX.OBS_INDEX_E{ind_obs,2} = ind_obs;
			
		end	
		
		
	end
	
	if strcmp(header(i-1,61:79),   'SYS / # / OBS TYPES') % Case of RINEX V3.0
		if (~strcmp(header(i-1,1), ' ') && ~strcmp(header(i-1,1),constellation)) % non-white char and other constellation
			constellation = header(i-1,1); % new constellation
			index_temp = 1;
			tab_obs = strsplit(strtrim(header(i-1,8:60)),' ');
			
			if strcmp(constellation,'G')
				%filling of RNX_OBS index
				for j=1:size(tab_obs,2)
					ind_obs = size(RNX.OBS_INDEX_G,1)+1;
					RNX.OBS_INDEX_G{ind_obs,1} = tab_obs{1,j};
					RNX.OBS_INDEX_G{ind_obs,2} = index_temp;
					index_temp = index_temp + 1;
				end
			elseif strcmp(constellation,'R')
				%filling of RNX_OBS index
				for j=1:size(tab_obs,2)
					ind_obs = size(RNX.OBS_INDEX_R,1)+1;
					RNX.OBS_INDEX_R{ind_obs,1} = tab_obs{1,j};
					RNX.OBS_INDEX_R{ind_obs,2} = index_temp;
					index_temp = index_temp + 1;
				end
			elseif strcmp(constellation,'E')
				%filling of RNX_OBS index
				for j=1:size(tab_obs,2)
					ind_obs = size(RNX.OBS_INDEX_E,1)+1;
					RNX.OBS_INDEX_E{ind_obs,1} = tab_obs{1,j};
					RNX.OBS_INDEX_E{ind_obs,2} = index_temp;
					index_temp = index_temp + 1;
				end
			end
		else
			tab_obs = strsplit(strtrim(header(i-1,8:60)),' ');
			if strcmp(constellation,'G')
				%filling of RNX_OBS index
				for j=1:size(tab_obs,2)
					ind_obs = size(RNX.OBS_INDEX_G,1)+1;
					RNX.OBS_INDEX_G{ind_obs,1} = tab_obs{1,j};
					RNX.OBS_INDEX_G{ind_obs,2} = index_temp;
					index_temp = index_temp + 1;
				end
			elseif strcmp(constellation,'R')
				%filling of RNX_OBS index
				for j=1:size(tab_obs,2)
					ind_obs = size(RNX.OBS_INDEX_R,1)+1;
					RNX.OBS_INDEX_R{ind_obs,1} = tab_obs{1,j};
					RNX.OBS_INDEX_R{ind_obs,2} = index_temp;
					index_temp = index_temp + 1;
				end
			elseif strcmp(constellation,'E')
				%filling of RNX_OBS index
				for j=1:size(tab_obs,2)
					ind_obs = size(RNX.OBS_INDEX_E,1)+1;
					RNX.OBS_INDEX_E{ind_obs,1} = tab_obs{1,j};
					RNX.OBS_INDEX_E{ind_obs,2} = index_temp;
					index_temp = index_temp + 1;
				end
			end
		end
	end
	   
	if(strcmp(header(i-1,61:73),   'END OF HEADER'))
		break;
	end

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [dataG,dataR,dataE,Nepoch]=read_data_3(frnx,RNX,epoch_max)
%% function [dataG,dataR,dataE,Nepoch]=read_data_3(frnx,RNX,epoch_max) : RINEX data loading (RINEX v3.__)
%% 
%% Jacques Beilin - ENSG/DPTS - 2010-02-25
%% Clement Fontaine - 2013-10-14
%%
%% Input :
%% - frnx : rinex file handle (already openned)
%% - RNX : Rinex header structure
%% - epoch_max : number of epochs to be loaded (if -1 : load all data)
%%
%% Output :
%% - dataG,dataR,dataE : GPS, Glonass and Galileo data structures (3D matrix  32 satellites x n observables x 3600 epochs)
%%   observabls order [ mjd obs ], where obs is defined in RNX (header)
%% - Nepoch : number of epochs loaded
%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Allocation (3600 epochs max = 24h - 30s or 1h -1s) - can work with more data but is slower
max_size = 3600;

% GPS
nb_obsG = size(RNX.OBS_INDEX_G,1); % number of observable
dataG=zeros(32,nb_obsG+1,max_size); % mjd + observable

%Glonass
nb_obsR = size(RNX.OBS_INDEX_R,1); % number of observable
dataR=zeros(32,nb_obsR+1,max_size); % mjd + observable

%Galileo
nb_obsE = size(RNX.OBS_INDEX_E,1); % number of observable
dataE=zeros(32,nb_obsE+1,max_size); % mjd + observable
		
% data reading
mjd = 0;
Nepoch = 0;
	
while ~ feof(frnx)
	line = fgetl(frnx);
	
	% new epoch
	if (strcmp(line(1:1),'>') && (strcmp(line(32:32),'0') || strcmp(line(32:32),'1'))) 
		
		Nepoch = Nepoch + 1;
				
		% date reading
		[vec] = sscanf(line(1:35),'%s %d %d %d %d %d %f %d %d');
		y = vec(2);
		m = vec(3);
		d = vec(4);
		hh = vec(5);
		mm = vec(6);
		ss = vec(7);
		[tgps]=ymdhms_t(y,m,d,hh,mm,ss);
		mjd=tgps.mjd;
		
		nb_sat = vec(9); % number of visible satellites at this epoch
		
		% loop on satellites
		for i = 1:nb_sat
				   
			line = fgetl(frnx);
						
			constellation = line(1:1);
			id = str2num(line(2:3));
			
			line = line(4:end);
			
			% Size of the line depends on the constellation			
			
			% GPS satellite 
			if strcmp(constellation,'G')
			
				% creates one observation line by satellite            
				str_obs = repmat(' ',1,16*nb_obsG);
				if(length(line)>0)
					str_obs(1:length(line)) = line;
				end
					
				% filling Vobs	
				[Vobs]=get_obs_from_str(str_obs,nb_obsG);
				Vobs=[mjd,Vobs];
					
				dataG(id,:,Nepoch) = Vobs;	
				
			% GLONASS satellite
			elseif strcmp(constellation,'R')
			
				% creates one observation line by satellite            
				str_obs = repmat(' ',1,16*nb_obsR);
				if(length(line)>0)
					str_obs(1:length(line)) = line;
				end
			
				% filling Vobs	
				[Vobs]=get_obs_from_str(str_obs,nb_obsR);
				Vobs=[mjd,Vobs];
					
				dataR(id,:,Nepoch) = Vobs;	

			% Galileo satellite
			elseif strcmp(constellation,'E')
			
				% creates one observation line by satellite            
				str_obs = repmat(' ',1,16*nb_obsE);
				if(length(line)>0)
					str_obs(1:length(line)) = line;
				end
				
				% filling Vobs	
				[Vobs]=get_obs_from_str(str_obs,nb_obsE);
				Vobs=[mjd,Vobs];
					
				dataE(id,:,Nepoch) = Vobs;	
			
			end
		         
		end
				
	end
	
	if(Nepoch>=epoch_max && epoch_max>-1) % breaks if epoch_max epochs have been read
		break;
	end

end

% resize matrix taking into account the number on epochs
dataG(:,:,Nepoch+1:end)=[];
dataR(:,:,Nepoch+1:end)=[];
dataE(:,:,Nepoch+1:end)=[];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [dataG,dataR,dataE,Nepoch]=read_data_2_11(frnx,RNX,epoch_max) %
%% function [dataG,dataR,dataE,Nepoch]=read_data_2_11(frnx,RNX,epoch_max) : RINEX data loading (RINEX v2.11)
%% 
%% Jacques Beilin - ENSG/DPTS - 2010-02-25
%% Clement Fontaine - 2013-10-14
%%
%% Input :
%% - frnx : rinex file handle (already openned)
%% - RNX : Rinex header structure
%% - epoch_max : number of epochs to be loaded (if -1 : load all data)
%%
%% Output :
%% - dataG,dataR,dataE : GPS, Glonass and Galileo data structures (3D matrix  32 satellites x n observables x 3600 epochs)
%%   observabls order [ mjd obs ], where obs is defined in RNX (header)
%% - Nepoch : number of epochs loaded
%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Allocation (3600 epochs max = 24h - 30s or 1h -1s) - works with more data but is slower
max_size = 3600;

nb_obs = size(RNX.OBS_INDEX_G,1);
nline_obs = 1 + floor(nb_obs / 5);

% GPS
dataG=zeros(32,nb_obs+1,max_size);

%Glonass
dataR=zeros(32,nb_obs+1,max_size);

%Galileo
dataE=zeros(32,nb_obs+1,max_size);
	
% nb of epochs	
Nepoch=0;

% data reading
mjd = 0;


while ~ feof(frnx)
	% epoch header reading
	line = fgetl(frnx);
	epoch_flag = sscanf(line(29:29),'%d');
	if (epoch_flag==0 || epoch_flag==1)
		Nepoch=Nepoch+1;
		
		% date reading
		[vec] = sscanf(line(1:26),'%d %d %d %d %d %f');
		y = vec(1);
		m = vec(2);
		d = vec(3);
		hh = vec(4);
		mm = vec(5);
		ss = vec(6);
		[tgps]=ymdhms_t(y,m,d,hh,mm,ss);
		mjd=tgps.mjd;
					
	    % searchs satellite list
		Nsat = sscanf(line(30:32),'%d');
		str_liste_sat=regexprep(line(33:length(line)),'(\s+$)','');
		Nline_header_epoch = ceil(Nsat / 12);
		if (Nline_header_epoch>1)
			% if Nsat>12, next line
			for i=2:Nline_header_epoch 
				line = fgetl(frnx);
				str_liste_sat=strcat(str_liste_sat,regexprep(line(33:length(line)),'(\s+$)',''));
			end
		end
		    
		% satellite data reading
		for i=1:Nsat
				
			%~ % creates one observation line by satellite            
			str_obs = repmat(' ',1,16*nb_obs);
			
			for j=1:nline_obs
				
				str_new = fgetl(frnx);
				len_str = length(str_new);
				
				if(len_str>0)
					str_obs(80*(j-1)+1:80*(j-1)+1 +length(str_new)-1) = str_new;
				end
		
			end

			%printf('>%s<\n',str_obs);
		
			[Vobs]=get_obs_from_str(str_obs,nb_obs);
			type_sat = str_liste_sat((i-1)*3+1:(i-1)*3+1);
			[num_sat] = sscanf(str_liste_sat((i-1)*3+2:(i-1)*3+3),'%d');
			
			% Filling obs
			
			% GPS satellite
			if (type_sat=='G')
				dataG(num_sat,:,Nepoch)=[mjd Vobs];

			% Glonass satellite
			elseif (type_sat=='R')
				dataR(num_sat,:,Nepoch)=[mjd Vobs]; 
				
			% Galileo satellite
			elseif (type_sat=='E')
				dataE(num_sat,:,Nepoch)=[mjd Vobs]; 
			
			end
								
		end
					
	else % suppress invalid epochs (epoch flag in [2:5])
		%%toc
		Ncomment=sscanf(line(30:32),'%d'); % number of lines to skip
		for i=1:Ncomment
			line = fgetl(frnx);
		end
	end
	
	
	if(Nepoch>=epoch_max && epoch_max>-1) % breaks if epoch_max epochs have been read
		break;
	end


end

% resize matrix taking into account the number on epochs
dataG(:,:,Nepoch+1:end)=[];
dataR(:,:,Nepoch+1:end)=[];
dataE(:,:,Nepoch+1:end)=[];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [Vobs]=get_obs_from_str(str,nb_obs)
%% function [Vobs]=get_obs_from_str(str,nb_obs)
%% Get obs from a line 
%%
%% Jacques Beilin - ENSG/DPTS - 2010-02-25
%% Clement Fontaine - 2013-10-14
%%
%% Input :
%% - str : string containing obs
%% - nb_obs : number of observations to be read
%%
%% Output :
%% - Vobs : observation vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vobs=zeros(1,nb_obs);

for j=1:nb_obs
		[x,c] = sscanf(str((j-1)*16+1:(j-1)*16+14),'%f');
		if c>0  ; Vobs(j) = x ; end                 
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [RNX_data] = qc(RNX_data)
%% function [RNX_data] = qc(RNX_data)
%% Data quality check
%%
%% Clement Fontaine - 2013-11-27
%%
%% Input :
%% - RNX_data : RNX_data structure
%%
%% Output :
%% - RNX_data : cleaned RNX_data structure
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RNX_dataG = RNX_data.G;
RNX_dataR = RNX_data.R;
RNX_dataE = RNX_data.E;

% Redim matrix -> pb if nb_epoch > 3600
nb_ep_max = max([size(RNX_dataG,3);size(RNX_dataR,3);size(RNX_dataE,3)]);

RNX_dataG = cat(3,RNX_dataG,zeros(size(RNX_dataG,1),size(RNX_dataG,2),nb_ep_max - size(RNX_dataG,3)));
RNX_dataR = cat(3,RNX_dataR,zeros(size(RNX_dataR,1),size(RNX_dataR,2),nb_ep_max - size(RNX_dataR,3)));
RNX_dataE = cat(3,RNX_dataE,zeros(size(RNX_dataE,1),size(RNX_dataE,2),nb_ep_max - size(RNX_dataE,3)));

RNX_data.G = RNX_dataG;
RNX_data.R = RNX_dataR;
RNX_data.E = RNX_dataE;

%% TODO

%~ % GPS
%~ for i = 1:32
%~ 
	%~ sat = squeeze(RNX_data.G(i,:,:));
	%~ if(sum(sum(sat(2:end,:)))==0)
		%~ RNX_data.G(i,:,:) = zeros(size(RNX_data.G,2),size(RNX_data.G,3));
	%~ end
%~ 
%~ end
%~ 
%~ % Glonass
%~ for i = 1:32
%~ 
	%~ sat = squeeze(RNX_data.R(i,:,:));
	%~ if(sum(sum(sat(2:end,:)))==0)
		%~ RNX_data.R(i,:,:) = zeros(size(RNX_data.R,2),size(RNX_data.R,3));
	%~ end
%~ 
%~ end
%~ 
%~ % Galileo
%~ for i = 1:32
%~ 
	%~ sat = squeeze(RNX_data.E(i,:,:));
	%~ if(sum(sum(sat(2:end,:)))==0)
		%~ RNX_data.E(i,:,:) = zeros(size(RNX_data.E,2),size(RNX_data.E,3));
	%~ end
%~ end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

