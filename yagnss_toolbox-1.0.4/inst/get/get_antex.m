function [ATX] = get_antex(ATX_header, ATX_data, ATX_name, ATX_freq, mjd)
%% function [ATX] = get_antex(ATX_header, ATX_data, ATX_name, ATX_freq, mjd)
%% Returns ATX structure corresponding to ATX_name antenna or satellite
%%
%% Clement Fontaine 2013-12-27
%%
%% Input : 
%% - ATX_header : structure containing header of ATX file (loaded by load_antex())
%% - ATX_data : structure containing data of ATX file (loaded by load_antex())
%% - ATX_name : name of antenna
%% - ATX_freq : freq used ('G01', 'G02', 'R01', 'R02')
%% - mjd : Modified Julian Date (mandatory if ATX_name corresponds to a satellite, otherwise useless) 
%%
%% Output : 
%% - ATX : structure containing antenna (or satellite) calibration 
%%
%%      TYPE = SVN                     : SVN = satellite, REC = receiver
%%      NAME = G03                     : antenna name
%%      FREQ = G01                     : freq in input
%%      CONST = G                      : constellation
%%      FREQ_CODE =  1                 : freq code
%%      NORTH =  0.279000000000000     : North offset (for receiver), or X Direction eccentricity of the satellite antenna phase center, relative to the satellite center of mass(m)
%%      EAST = 0                       : East offset (for receiver), or Y Direction eccentricity of the satellite antenna phase center, relative to the satellite center of mass (m)
%%      UP =  2.79260000000000         : Up offset (for receiver), or Z Direction eccentricity of the satellite antenna phase center, relative to the satellite center of mass (m)
%%      NOAZI =                        : NOn-AZImuth dependant phase pattern (m)
%%
%%        -0.00080000000000000
%%        -0.00090000000000000
%%        -0.00090000000000000
%%        -0.00080000000000000
%%        -0.00040000000000000
%%         0.00020000000000000
%%         0.00080000000000000
%%         0.00130000000000000
%%         0.00140000000000000
%%         0.00120000000000000
%%         0.00070000000000000
%%         0.00000000000000000
%%        -0.00040000000000000
%%        -0.00070000000000000
%%        -0.00090000000000000
%%        -0.00090000000000000
%%        -0.00090000000000000
%%        -0.00090000000000000
%%
%%      ETAL_AZI = [](0x0)             : Azi phase pattern (for receiver only) (m)
%%      VAZI = [](0x0)                 : Azimuthal angles (deg) (for receiver only)
%%      DAZI =                         : Zenithal angles (deg)
%%
%%          0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% output
ATX = struct;
ATX.TYPE = '';
ATX.NAME = '';
ATX.FREQ = '';

ATX.CONST = '';
ATX.FREQ_CODE = 0;
ATX.NORTH = 0;
ATX.EAST = 0;
ATX.UP = 0;
ATX.NOAZI = [];
ATX.ETAL_AZI = [];
ATX.VAZI = [];
ATX.DAZI = [];

if ~isfield(ATX_header,'NAME_INDEX')
	tool_print_info('Structure ATX is empty',3);
	return;
end


% looking for antenna

% Satellite
if(nargin==5)

	Nant = find(ismember(ATX_header.NAME_INDEX,ATX_name)==1);
	Nant2 = [];
	
	for i = 1:length(Nant)
		
		if strcmp(ATX_data{Nant(i),1}.TYPE,'SVN') % is a sat

			if isfield(ATX_data{Nant(i),1},'UNTIL') % not last generation
			
				if(mjd>=ATX_data{Nant(i),1}.FROM && mjd<=ATX_data{Nant(i),1}.UNTIL)
					
					Nant2 = [Nant2;Nant(i)];
					
				end
			
			else % last generation
			
				if(mjd>=ATX_data{Nant(i),1}.FROM)
				
					Nant2 = [Nant2;Nant(i)];
					
				end
			
			end

		else
		
			continue;
			
		end
		
	end
	
	Nant = Nant2;

% Receiver antenna
else

	Nant = find(ismember(ATX_header.NAME_INDEX,ATX_name)==1);

end

if length(Nant) == 0
	tool_print_info(sprintf('Unable to find antenna %s !\n',ATX_name),3)
	return;
end


% Frequency

Nfreq = -1;
for i=2:ATX_data{Nant,1}.N_FREQUENCIES+1
	[S, E, TE, M, T, NM] = regexp (ATX_data{Nant,i}.FREQUENCY_NAME, ATX_freq);
	if S>=0 
		Nfreq = i;	
	end	
end
if Nfreq <0
	tool_print_info(sprintf('Unable to find frequency %s for antenna %s !\n',ATX_freq,ATX_name),3)
	return;
end


ATX.TYPE = ATX_data{Nant,1}.TYPE;
ATX.NAME = ATX_data{Nant,1}.NAME;
ATX.FREQ = ATX_data{Nant,Nfreq}.FREQUENCY_NAME;

ATX.CONST = ATX_data{Nant,Nfreq}.SAT_SYSTEM;
ATX.FREQ_CODE = ATX_data{Nant,Nfreq}.FREQ_CODE;

ATX.NORTH = ATX_data{Nant,Nfreq}.NORTH*1e-3;
ATX.EAST = ATX_data{Nant,Nfreq}.EAST*1e-3;
ATX.UP = ATX_data{Nant,Nfreq}.UP*1e-3;

ATX.NOAZI = ATX_data{Nant,Nfreq}.NOAZI*1e-3;
ATX.DAZI = ATX_data{Nant,1}.ZEN1:ATX_data{Nant,1}.DZEN:ATX_data{Nant,1}.ZEN2;

if nargin ~= 5
	ATX.ETAL_AZI = ATX_data{Nant,Nfreq}.ETAL_AZI*1e-3;
	ATX.VAZI = ATX_data{Nant,Nfreq}.VAZI;
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
