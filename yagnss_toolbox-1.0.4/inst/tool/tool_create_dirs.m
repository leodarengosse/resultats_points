function [err] = tool_create_dirs(dirs)
%% function [err] = tool_create_dirs(dirs)
%%
%% Create directories
%%
%% Clement Fontaine 2013-11-18
%%
%% Input : 
%% - dirs : directory to create
%%		ex : './ee/aa/pp'
%%
%% Output : 
%% - err : error code (now set to 0 ... TODO)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Windows or Linux separator
sep = filesep();

if ~strcmp(dirs,'')

	dir = strsplit(dirs,sep);
	if strcmp(dir(1),'.')
		dir(1) = [];
	end
		
	curr_dir = dir{1};
	mkdir(curr_dir);

	for i = 1:length(dir)-1

		if mkdir(curr_dir,dir{i+1})		
			tool_print_info(sprintf('Directory : %s created in %s directory',dir{i+1},curr_dir()),1);
		end
		
		curr_dir = strcat(curr_dir,sep,dir{i+1});
				
	end

end

err=0;

end
