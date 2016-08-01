function [ATX_header,ATX_data]=load_antex(filename)
%% function [ATX_header,ATX_data]=load_antex(filename)
%% Loading ANTEX file
%% 
%% Jacques Beilin - ENSG/DPTS - 2011-12-26
%%
%% Input :
%% - filename : ANTEX file
%% Output :
%% - ATX_header : structure containing ANTEX header
%% - ATX_data : cell tab conaining ANTEX data
%%
%% Use get_antex() to extract the right atx
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bin_filename = sprintf('%s.bin',filename);

[Fid,msg]=fopen(bin_filename);
%~ Fid = -4;
if Fid>-1
%if Fid==NaN

    tool_print_info(sprintf ('Loading binary file  %s (Original ANTEX file %s)',bin_filename,filename),1);
    fclose(Fid); 
    
    % load
    if exist('OCTAVE_VERSION') % octave 
		load ('-binary', bin_filename, 'ATX_header','ATX_data');
	else % matlab
		load (bin_filename, '-mat', 'ATX_header','ATX_data');
	end
	

else

    % initialization
    ATX_header=[];
    ATX_data=[];
    ATX_header.VERSION=0;    

    % File opening
    [f,mess] = fopen(filename, 'rt');

    % check if file exists
	if (f<=-1)
		tool_print_info(sprintf('Unable to open %s',filename),3);
		return;
	else

        % header reading
        [ATX_header]=read_header_atx(f,ATX_header);

        % data reading
        [ATX_data]=read_atx_data(f,ATX_header);
        
        % set index
        [ATX_header] = set_ATX_NAME_INDEX(ATX_header, ATX_data);

        status = fclose(f);

		% save
		if exist('OCTAVE_VERSION') % octave 
			save('-binary', bin_filename, 'ATX_header','ATX_data');
		else % matlab
			save(bin_filename,'-mat', 'ATX_header','ATX_data');      
		end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ATX]=read_header_atx(f,ATX);
%% function [ATX_header]=read_header_atx(f,ATX_header)
%% 
%% Jacques Beilin - ENSG/DPTS - 2011-12-26
%%
%% Input :
%% - f : file_id (already opened)
%% - ATX : structure containing ATX (empty)
%%
%% Output :
%% - ATX : structure containing ATX
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ATX.VERSION=0;
    i=1;
    while(i<100) 
        line=fgetl(f);
        header(i,1:length(line)) = line;
        
        if(strcmp(header(i,61:73),   'END OF HEADER'))
			break;
		end
		
        i = i + 1;
        if strcmp(header(i-1,61:80),   'ANTEX VERSION / SYST') 
            ATX.VERSION=str2double(header(i-1,1:10));
            ATX.TYPE=header(i-1,21:21) ; 
        end
    end 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [ATX]=read_atx_data(f,ATX_header);
%% function [ATX_data]=read_atx_data(f,ATX_header) 
%% ATX data reading
%%
%% Jacques Beilin - ENSG/DPTS - 2011-12-26
%%
%% Input :
%% - f :  file_id (already opened)
%% - ATX_header : structure containing ATX header
%% Output :
%% - ATX : structure containing ATX
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% allocation for 2880 epochs (24h - 30s)
default_NAntenna = 4;
ATX = cell(default_NAntenna,10);
% data reading
Nant=0;
while ~ feof(f)
    line = fgetl(f);
    if strcmp(line(61:76),   'START OF ANTENNA')
		Nant=Nant+1;
		[ATX]=read_antenna(f,ATX,Nant);
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [ATX]=read_antenna(f,ATX,Nant);
%% function [ATX_data]=read_antenna(f,ATX_header)
%% Antenna data reading
%% 
%% Jacques Beilin - ENSG/DPTS - 2011-12-26
%%
%% Input :
%% - f : file_id (already opened)
%% - ATX_header : structure containing ATX header
%% Output :
%% - ATX_data : structure containing ATX data
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% allocation for 2880 epochs (24h - 30s)
TYPE = '';
while ~ feof(f)
    line = fgetl(f);
    
    if strcmp(line(61:74),'END OF ANTENNA')
        %ATX{Nant,[1:4]}
		return;
	end	
    
	if strcmp(line(61:76),'TYPE / SERIAL NO')
		if regexp(line(21:21),'[GRE]')==1
			sat = 1;
			ATX{Nant,1}.TYPE='SVN';
		else 
			sat = 0;
			ATX{Nant,1}.TYPE='REC';
		end
		
		if sat==1 
			ATX{Nant,1}.GEN = line(1:20);
			ATX{Nant,1}.PRN = sscanf(line(21:40),'%s');
			ATX{Nant,1}.SVN = sscanf(line(41:50),'%s');
			SN =  sscanf(line(51:60),'%s');
			ATX{Nant,1}.NAME=ATX{Nant,1}.PRN;
		else
			ATX{Nant,1}.SN =  sscanf(line(21:40),'%s');
			ATX{Nant,1}.NAME=line(1:20);

		end
	end
	
	if strcmp(line(61:64),'DAZI')
		ATX{Nant,1}.DAZI=sscanf(line(3:8),'%f');
	end
	
	if strcmp(line(61:78),'ZEN1 / ZEN2 / DZEN')
		ATX{Nant,1}.ZEN1=sscanf(line(3:8),'%f');
		ATX{Nant,1}.ZEN2=sscanf(line(9:14),'%f');
		ATX{Nant,1}.DZEN=sscanf(line(15:20),'%f');
	end
	
	if strcmp(line(61:76),'# OF FREQUENCIES')
		ATX{Nant,1}.N_FREQUENCIES=sscanf(line(1:6),'%d');
		ATX{Nant,1}.N_FREQ=ATX{Nant,1}.N_FREQUENCIES;
	end
	
	if strcmp(line(61:70),'VALID FROM')
		[vec] = sscanf(line(1:60),'%d %d %d %d %d %f');
		y = vec(1);
		m = vec(2);
		d = vec(3);
		hh = vec(4);
		mm = vec(5);
		ss = vec(6);
		[tgps]=ymdhms_t(y,m,d,hh,mm,ss);
		mjd=tgps.mjd;
		ATX{Nant,1}.FROM=mjd;
	end
	
	if strcmp(line(61:71),'VALID UNTIL')
		[vec] = sscanf(line(1:60),'%d %d %d %d %d %f');
		y = vec(1);
		m = vec(2);
		d = vec(3);
		hh = vec(4);
		mm = vec(5);
		ss = vec(6);
		[tgps]=ymdhms_t(y,m,d,hh,mm,ss);
		mjd=tgps.mjd;
		ATX{Nant,1}.UNTIL=mjd;
	end
	
	if strcmp(line(61:70),'SINEX CODE')
		ATX{Nant,1}.SINEX_CODE=sscanf(line(1:10),'%s');
	end
	 
	if strcmp(line(61:78),'START OF FREQUENCY')
		Nfreq = ATX{Nant,1}.N_FREQUENCIES - ATX{Nant,1}.N_FREQ + 2;
		ATX{Nant,Nfreq}.FREQUENCY_NAME=sscanf(line(4:6),'%s');
		ATX{Nant,Nfreq}.SAT_SYSTEM=sscanf(line(4:4),'%s');
		ATX{Nant,Nfreq}.FREQ_CODE=sscanf(line(5:6),'%d');
		[ATX]=read_frequency(f,ATX,Nant,Nfreq);
		ATX{Nant,1}.N_FREQ = ATX{Nant,1}.N_FREQ - 1;	
	end

	
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ATX]=read_frequency(f,ATX,Nant,Nfreq);
%% function [ATX]=read_frequency(f,ATX,Nant,Nfreq); 
%% 
%% 
%% Jacques Beilin - ENSG/DPTS - 2011-12-26
%%
%% Input :
%% - f : file_id (already opened)
%% - ATX : structure containing ATX
%% - Nant : current antenna number
%% - Nfreq : calibrated frequecy number
%% Output :
%% - ATX : structure containing ATX (with a new frequency)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TYPE = '';
while ~ feof(f)
    line = fgetl(f);
    
    if strcmp(line(61:76),'END OF FREQUENCY')
		return
	end	
	
	if strcmp(line(61:77),'NORTH / EAST / UP')
		ATX{Nant,Nfreq}.NORTH=sscanf(line(1:10),'%f');
		ATX{Nant,Nfreq}.EAST=sscanf(line(11:20),'%f');
		ATX{Nant,Nfreq}.UP=sscanf(line(21:30),'%f');
	end
	
	if strcmp(line(4:8),'NOAZI')
		m = (ATX{Nant,1}.ZEN2 - ATX{Nant,1}.ZEN1) / ATX{Nant,1}.DZEN + 1; 
		vec=sscanf(line([9:end]), '%f', m);
		ATX{Nant,Nfreq}.NOAZI=vec;
		
		if ATX{Nant,1}.DAZI>0
			nlignes_azi = 360 / ATX{Nant,1}.DAZI+1;
			ATX{Nant,Nfreq}.ETAL_AZI = zeros(nlignes_azi,m);
			
			for i=1:nlignes_azi
			    line = fgetl(f);
				ATX{Nant,Nfreq}.VAZI(i)=sscanf(line(1:8),'%f');
				vec=sscanf(line([9:end]), '%f', m);
				ATX{Nant,Nfreq}.ETAL_AZI(i,:)=vec';
			end
				
		end
	end
	 
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ATX_header] = set_ATX_NAME_INDEX(ATX_header, ATX_data)
%% function [ATX_header] = set_ATX_NAME_INDEX(ATX_header, ATX_data)
%%
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NAME_INDEX = cell(size(ATX_data,1),1);
for i = 1:size(ATX_data,1)
	NAME_INDEX{i} = ATX_data{i,1}.NAME;
end

ATX_header.NAME_INDEX = NAME_INDEX;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
