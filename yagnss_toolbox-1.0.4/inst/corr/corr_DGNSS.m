function [PR_correction] = corr_DGNSS(PRC,const,PRN)
%% function [PR_correction] = corr_DGNSS(PRC,const,PRN)
%% Calculate a DGNSS correction from PRC structure
%%
%% Clement Fontaine - 2013-11-13
%%
%% Input
%% - PRC : structure :
%%		- corr : PRC correction matrix 
%%		- index : sat_index (cell tab of sat id (format 'A1I2' : ex 'G01'))    
%%
%% - const : constellation 
%% - PRN : satellite id
%%
%% Output
%% - PR_correction : DGNSS correction for given satellite at a given date (m)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Output
PR_correction = 0;

sat = sprintf('%s%02d',const,PRN);
sat_pos = find(ismember(PRC.index,sat));

if sat_pos>0
	PR_correction = PRC.corr(sat_pos);
else
	tool_print_info(sprintf('DGNSS does not contain PRC of %s: PRC set to 0',sat),2);
	return;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
