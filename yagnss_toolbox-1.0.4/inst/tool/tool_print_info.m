function tool_print_info(str,code)
%% function tool_print_info(str,code)
%%
%% Print strings in terminal and in a log file
%%
%% Clement Fontaine 2013-11-18
%%
%% Input : 
%% - str : string to print
%% - code : code
%%		- 0 : print in screen, not in file
%%		- 1 : normal
%%		- 2 : warning
%%		- 3 : error
%%		- 4 : print in file, not in screen
%%
%% Global variables used : 
%% - LOG_FILE : file name of log. If not defined, no log file
%% - VERBOSE : id VERBOSE == 0 : no print in terminal
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global LOG_FILE;
global VERBOSE;

if size(LOG_FILE)>0

	file_id = fopen(LOG_FILE,'a');
	
	if file_id >-1
		
		switch code
		
			% normal
			case 1
				fprintf(file_id,'%s\n',str);
			
			% warning
			case 2
				fprintf(file_id,'WARNING : %s\n',str);
			
			% error
			case 3
				fprintf(file_id,'ERROR : %s\n',str);
			
			% print in file, not in screen
			case 4
				fprintf(file_id,'%s\n',str);
	
		end
		
		fclose(file_id);
	
	end

end

if ~VERBOSE == 0

	switch code
	
		% print in screen, not in file
		case 0
			fprintf('%s\n',str);
	
		% normal
		case 1
			fprintf('%s\n',str);		
			
		% warning
		case 2	
			fprintf('WARNING : %s\n',str);	
			
		% error
		case 3
			fprintf('ERROR : %s\n',str);
	end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
