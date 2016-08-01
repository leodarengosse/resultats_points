function [NAV_header,NAV_data]=load_rinex_n(filenames)
%% function [NAV_header,NAV_data]=load_rinex_n(filenames)
%% navigation RINEX loading
%%
%% Jacques Beilin - ENSG/DPTS - 2010-11-09
%% Clement Fontaine - 2013-10-15
%%
%% Input :
%% - filenames : navigation files in a cell array (RINEX v2.11 or RINEX v3.__)
%%
%% Output :
%% - NAV_header : structure contaning header
%%   Contents : parameters of Klobuchar correction, corrections to transform
%%              the system time to UTC or other time systems
%%   NAV_header =
%%   {
%%   	VERSION =  3.02000000000000
%%      TYPE = N
%%      LEAP_SECONDS =  16
%%		GPSA = 0   0   0   0    % ION ALPHA for RINEX V2.11
%%		GPSB = 0   0   0   0    % ION BETA for RINEX V2.11
%%		GAL = 0   0   0   0     % Ionospheric correction for Galileo
%%      GPUT = % GPS to UTC
%%      etc ... (cf doc RINEX 3.02)
%%
%%   }
%% - NAV_data : structure containing 3D matrix containing navigation data for GPS, GLONASS and GALILEO
%%   NAV_data.G, NAV_data.R, NAV_data.E
%%   Index order :
%%   	- Satellite number (PRN), between 1 and 32
%%   	- Navigation message data (first is mjd)
%%   	- Epoch (Warning : matrix may contain several times the same data)
%%
%%   INFO  :
%%   Use the function get_ephemeride to get the data of navigation files
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Separator for windows or linux
sep = filesep();

if ischar(filenames)
    filenames = cellstr(filenames);
end


% bin filename

% directory
str = strsplit(filenames{1}(1:end-4),sep);
ext = filenames{1}(end-2:end-1);
dir = '';

for j = 1: length(str)-1
    dir = sprintf('%s%s',dir,str{j},sep);
end

bin_filenames = sprintf('%s%s.%sp.bin',dir,str{end},ext);



[Fid,msg]=fopen(bin_filenames);

% debug
%Fid = -1; % do not read .bin even if it exists
% end debug

if Fid>-1
    %if Fid==NaN
    
    filenames_temp = '';
    for i = 1:length(filenames)
        filenames_temp = sprintf('%s %s',filenames_temp,filenames{i});
    end
    
    tool_print_info(sprintf ('Loading binary file  %s (Original RINEX file %s)',bin_filenames,filenames_temp),1);
    fclose(Fid);
    
    % load
    if exist('OCTAVE_VERSION') % octave
        load('-binary', bin_filenames, 'NAV_header','NAV_data');
    else % matlab
        load (bin_filenames,'-mat', 'NAV_header','NAV_data');
    end
    
else
    
    % Initialization
    NAV_header=[];
    NAV_data=[];
    NAV_dataG=[];
    NAV_dataR=[];
    NAV_dataE=[];
    NAV_dataG_temp=zeros(1,1,1);
    NAV_dataR_temp=zeros(1,1,1);
    NAV_dataE_temp=zeros(1,1,1);
    
    NAV_header.VERSION=0;
    
    for n_fich = 1:length(filenames)
        
        % File opening
        [frnx,mess] = fopen(filenames{n_fich}, 'rt');
        
        % check if file exists
        if (frnx<=-1)
            tool_print_info(sprintf('Unable to open %s',filenames{n_fich}),3);
            return;
        else
            
            tool_print_info(sprintf ('Loading RINEX navigation file %s',filenames{n_fich}),1);
            
            % header loading
            [NAV_header]=read_header_Nav(frnx,NAV_header);
            
            % navigation message loading
            if NAV_header.VERSION>=3.0
                [NAV_dataG_temp,NAV_dataR_temp,NAV_dataE_temp]=read_nav_GPS_3(frnx);
            else %
                [NAV_dataG_temp]=read_nav_GPS_211(frnx);
            end
            
            status = fclose(frnx);
            
            % affectation
            
            if sum(sum(NAV_dataG_temp(:,1,:)))>0
                NAV_data.G=NAV_dataG_temp;
            end
            
            if sum(sum(NAV_dataR_temp(:,1,:)))>0
                NAV_data.R=NAV_dataR_temp;
            end
            
            if sum(sum(NAV_dataE_temp(:,1,:)))>0
                NAV_data.E=NAV_dataE_temp;
            end
            
            
        end
        
    end
    
    % save
    if exist('OCTAVE_VERSION') % octave
        save('-binary', bin_filenames, 'NAV_header','NAV_data');
    else % matlab
        save(bin_filenames,'-mat', 'NAV_header','NAV_data');
    end
    
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





function [RNX]=read_header_Nav(frnx,RNX)
%% function [RNX]=read_header(frnx,RNX)
%% Header loading
%%
%% Jacques Beilin - ENSG/DPTS - 2010-02-25
%% Clement Fontaine - 2013-10-14
%%
%% Input :
%% - frnx : RINEX file handle (already openned)
%% - RNX : NAV_header structure
%%
%% Output :
%% - RNX : RINEX header structure
%%
%% Fields 'TIME SYSTEM CORR' depend on the content of header : test if a field exists = isfield(RNX,'GAUT')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(RNX,'VERSION')
    RNX.VERSION=0;
end

if ~isfield(RNX,'TYPE')
    RNX.TYPE='';
end

if ~isfield(RNX,'LEAP_SECONDS')
    RNX.LEAP_SECONDS=0;
end

if ~isfield(RNX,'GPSA')
    RNX.GPSA=zeros(1,4); % = ION ALPHA (GPS)
end

if ~isfield(RNX,'GPSB')
    RNX.GPSB=zeros(1,4); % = ION BETA (GPS)
end

if ~isfield(RNX,'GAL')
    RNX.GAL=zeros(1,4); % Ionospheric correction for Galileo
end

% Depends on rinex file : test if a field exists = isfield(RNX,'GAUT')
%RNX.GAUT=zeros(1,4); % GAL to UTC
%RNX.GPUT=zeros(1,4); % GPS to UTC
%RNX.GLUT=zeros(1,4); % GLO to UTC
%RNX.GPGA=zeros(1,4); % GPS to Gal
%RNX.GPGL=zeros(1,4); % GLO to UTC

i=1;
while(i<100)
    line=fgetl(frnx);
    header(i,1:length(line)) = line;
    i=i+1;
    

    % RINEX version number and type
    if strcmp(header(i-1,61:80),   'RINEX VERSION / TYPE')
        RNX.VERSION=str2double(header(i-1,1:10));
        RNX.TYPE=header(i-1,21:21) ;
    end
    
    % Leap seconds
    if strcmp(header(i-1,61:72),   'LEAP SECONDS')
        RNX.LEAP_SECONDS=str2num(header(i-1,1:6));
    end

    if RNX.VERSION>=3
        
        % Ionospheric correction
        if strcmp(header(i-1,61:76),   'IONOSPHERIC CORR')
            
            if strcmp(header(i-1,1:4),   'GPSA')
                for j=1:4
                    RNX.GPSA(j)=str2num(header(i-1,6+(j-1)*12:6+(j)*12));
                end
            end
            
            if strcmp(header(i-1,1:4),   'GPSB')
                for j=1:4
                    RNX.GPSB(j)=str2num(header(i-1,6+(j-1)*12:6+(j)*12));
                end
            end
            
            if strcmp(header(i-1,1:3),   'GAL')
                for j=1:4
                    RNX.GAL(j)=str2num(header(i-1,6+(j-1)*12:6+(j)*12));
                end
            end
            
        end
        
        % Time system coord
        if strcmp(header(i-1,61:76),   'TIME SYSTEM CORR')
            
            % RNX.NOM_CHANGEMENT_TYPE(i) = val;
            TSC = header(i-1,1:4);
            RNX.(TSC) = [str2num(header(i-1,6:22)),str2num(header(i-1,23:38)),str2num(header(i-1,39:45)),str2num(header(i-1,46:50)),str2num(header(i-1,51:55))];
        end
        
    else % Rinex 2.x
    
        % Ionospheric correction
        if regexp(line,'ION\sALPHA')>0
            for j=1:4
            line(2+(j-1)*12:2+j*12)
                RNX.GPSA(j)=str2num(line(2+(j-1)*12:2+j*12));
            end
        end
        if regexp(line,'ION\sBETA')>0
            for j=1:4
                RNX.GPSB(j)=str2num(header(i-1,2+(j-1)*12:2+j*12));
            end
        end
        
        % Time system coord
        if strcmp(header(i-1,61:80),   'DELTA-UTC: A0,A1,T,W')
            RNX.GPUT(1)=str2num(header(i-1,3:22));
            RNX.GPUT(2)=str2num(header(i-1,23:41));
            RNX.GPUT(3)=str2num(header(i-1,42:50));
            RNX.GPUT(4)=str2num(header(i-1,51:59));
        end
    end
    
    if(strcmp(header(i-1,61:73),   'END OF HEADER'))
        break;
    end
    
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





function [NAV_data]=read_nav_GPS_211(frnx)
%% function [NAV_data]=read_nav_GPS_211(frnx,NAV_header) : Navigation message loading (RINEX v2.11)
%%
%% Jacques Beilin - ENSG/DPTS - 2010-02-25
%% Clement Fontaine - 2013-10-14
%%
%% Input :
%% - frnx : RINEX file handle (already openned)
%%
%% Output :
%% - NAV_data : navigation message structure
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NAV_data=zeros(32,28,24); % 32 satellites, 28 infos. Can work with more data but is slower

% data reading
mjd = 0;
Nepoch=zeros(32,1);
while ~ feof(frnx)
    % Satellite type and number (PRN)
    line = fgetl(frnx);
    PRN = sscanf(line(1:2),'%d');
    
    % date reading
    [vec] = sscanf(line(4:22),'%d %d %d %d %d %f');
    y = vec(1);
    m = vec(2);
    d = vec(3);
    hh = vec(4);
    mm = vec(5);
    ss = vec(6);
    [tgps]=ymdhms_t(y,m,d,hh,mm,ss);
    mjd=tgps.mjd;
    
    Nepoch(PRN)=Nepoch(PRN)+1;
    
    
    NAV_data(PRN,1,Nepoch(PRN))=mjd;
    
    index = 2;
    for i=1:3
        str=regexprep(line(23+(i-1)*19:22+i*19),'D','E');
        NAV_data(PRN,index,Nepoch(PRN))=sscanf(str,'%f');
        index=index+1;
    end
    
    for j=1:6
        line = fgetl(frnx);
        for i=1:4
            str=regexprep(line(4+(i-1)*19:3+i*19),'D','E');
            NAV_data(PRN,index,Nepoch(PRN))=sscanf(str,'%f');
            index=index+1;
        end
    end
    
    line = fgetl(frnx);
    str=regexprep(line(4:22),'D','E');
    NAV_data(PRN,index,Nepoch(PRN))=sscanf(str,'%f');
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





function [NAV_dataG,NAV_dataR,NAV_dataE]=read_nav_GPS_3(frnx)
%% function [NAV_data]=read_nav_GPS_3(frnx,NAV_header) : Navigation message loading (RINEX v3.__)
%%
%% Jacques Beilin - ENSG/DPTS - 2010-02-25
%% Clement Fontaine - 2013-10-14
%%
%% Input :
%% - frnx : RINEX file handle (already openned)
%%
%% Output :
%% - NAV_data : navigation message structure
%%


NAV_dataG=zeros(32,29,24); % GPS
NAV_dataR=zeros(32,16,24); % GLONASS
NAV_dataE=zeros(32,28,24); % GALILEO

mjd = 0;
NepochG=zeros(32,1);
NepochR=zeros(32,1);
NepochE=zeros(32,1);

jmaxG = 7;
jmaxR = 3;
jmaxE = 7;

imaxG=[4,4,4,4,4,4,1];
imaxR=[4,4,4];
imaxE=[4,4,4,4,3,4,1];

% data reading

while ~ feof(frnx)
    
    
    line = fgetl(frnx);
    % satellite id
    constellation = line(1:1);
    PRN = sscanf(line(2:3),'%d');
    
    if ~(strcmp(constellation,'G') || strcmp(constellation,'R') || strcmp(constellation,'E'))
        % loop until a new satellite
        continue;
    end
    
    % date reading
    [vec] = sscanf(line(5:23),'%d %d %d %d %d %f');
    y = vec(1);
    m = vec(2);
    d = vec(3);
    hh = vec(4);
    mm = vec(5);
    ss = vec(6);
    [tgps]=ymdhms_t(y,m,d,hh,mm,ss);
    mjd=tgps.mjd;
    
    % depends on the constellation
    
    % GPS
    if strcmp(constellation,'G')
        NepochG(PRN) = NepochG(PRN)+1;
        
        NAV_dataG(PRN,1,NepochG(PRN))=mjd;
        
        index = 2;
        for i=1:3
            
            str=regexprep(line(24+(i-1)*19:23+i*19),'D','E');
            NAV_dataG(PRN,index,NepochG(PRN))=sscanf(str,'%f');
            index=index+1;
            
        end
        
        for j=1:jmaxG
            line = fgetl(frnx);
            for i=1:imaxG(j)
                
                str=regexprep(line(5+(i-1)*19:4+i*19),'D','E');
                NAV_dataG(PRN,index,NepochG(PRN))=sscanf(str,'%f');
                index=index+1;
            end
        end
        
        
    elseif strcmp(constellation,'R')
        NepochR(PRN) = NepochR(PRN)+1;
        
        NAV_dataR(PRN,1,NepochR(PRN))=mjd;
        
        index = 2;
        for i=1:3
            
            str=regexprep(line(24+(i-1)*19:23+i*19),'D','E');
            NAV_dataR(PRN,index,NepochR(PRN))=sscanf(str,'%f');
            index=index+1;
            
        end
        
        for j=1:jmaxR
            line = fgetl(frnx);
            for i=1:imaxR(j)
                
                str=regexprep(line(5+(i-1)*19:4+i*19),'D','E');
                NAV_dataR(PRN,index,NepochR(PRN))=sscanf(str,'%f');
                
                index=index+1;
            end
        end
        
    elseif strcmp(constellation,'E')
        
        NepochE(PRN) = NepochE(PRN)+1;
        NAV_dataE(PRN,1,NepochE(PRN))=mjd;
        
        index = 2;
        for i=1:3
            
            str=regexprep(line(24+(i-1)*19:23+i*19),'D','E');
            NAV_dataE(PRN,index,NepochE(PRN))=sscanf(str,'%f');
            index=index+1;
            
        end
        
        for j=1:jmaxE
            line = fgetl(frnx);
            for i=1:imaxE(j)
                
                str=regexprep(line(5+(i-1)*19:4+i*19),'D','E');
                NAV_dataE(PRN,index,NepochE(PRN))=sscanf(str,'%f');
                
                index=index+1;
            end
        end
        
    end
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
