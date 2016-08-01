function [RNX_header1,RNX_data1,RNX_header2,RNX_data2] = tool_add_SA(RNX_header1,RNX_data1,RNX_header2,RNX_data2)
%% function tool_add_SA()
%% Add a virtual 'SA' in observations
%%
%% Input : 
%% - RNX_header1,RNX_data1,RNX_header2,RNX_data2 : rinex obs structure obtained from load_rinex_o()
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tps = 0;
const = 'GRE';
tic
% On parcours un RNX_data.G, et on cherche l'epoque correspondante dans l'autre
% rinex

% case 1 : nargin == 2 (one receiver)
if nargin == 2

	tool_print_info('',1);
	tool_print_info(sprintf('Add noise on %s',RNX_header1.MARKER_NAME),1);
	
	for const_i = 1:3
	
		obs_id1 = eval(strcat('RNX_header1.OBS_INDEX_',const(const_i)));
	
		for i = 1:32
		
			data = RNX_data1.(const(const_i));
			data = squeeze(data(i,:,:))';
			
			off = 300*randn(size(data,1),size(data,2)-1);
			data(:,2:end) = data(:,2:end)+off; 
			
			RNX_data1.(const(const_i))(i,:,:) = data';
		end
	end

% case 2 : nargin == 4 (two receivers)
else

	tool_print_info('',1);
	tool_print_info(sprintf('Add noise on %s and %s',RNX_header1.MARKER_NAME,RNX_header2.MARKER_NAME),1);

	for const_i = 1:3
	
		obs_id1 = eval(strcat('RNX_header1.OBS_INDEX_',const(const_i)));
		obs_id2 = eval(strcat('RNX_header2.OBS_INDEX_',const(const_i)));
		
		for i = 1:32
		
			data1 = RNX_data1.(const(const_i));
			data1 = squeeze(data1(i,:,:))';
			
			data2 = RNX_data2.(const(const_i));
			data2 = squeeze(data2(i,:,:))';
			
			off = 300*randn(size(data1,1),size(data1,2)-1);
			
			data1(:,2:end) = data1(:,2:end)+off; 
			data2(:,2:end) = data2(:,2:end)+off; 
			
			RNX_data1.(const(const_i))(i,:,:) = data1';
			RNX_data2.(const(const_i))(i,:,:) = data2';
		
		end
		
	end

end


%~ 
%~ 
%~ 
%~ 
%~ for const_i = 1:3
%~ 
	%~ const(const_i)
	%~ 
	%~ obs_id1 = eval(strcat('RNX_header1.OBS_INDEX_',const(const_i)));
	%~ obs_id2 = eval(strcat('RNX_header2.OBS_INDEX_',const(const_i)));
	%~ 
	%~ data_id = [];
	%~ for i = 1:size(obs_id1,1)
		%~ 
		%~ if (obs_id1{i,1}(1) == 'L' || obs_id1{i,1}(1) == 'P' || obs_id1{i,1}(1) == 'C')
		%~ if (obs_id1{i,1}(1) == 'P' || obs_id1{i,1}(1) == 'C') % only code
		%~ 
			%~ for j = 1:size(obs_id2,1)
			%~ 
				%~ if obs_id1{i,1} == obs_id2{j,1}
					%~ 
					%~ data_id = [data_id;i,j];
					%~ data_id = [data_id;i];
				%~ 
				%~ end
				%~ 
			%~ end
			%~ 
		%~ end
	%~ end
%~ 
	%~ for i = 1:32
		%~ i
		%~ 
		%~ 
		%~ data1 = RNX_data1.(const(const_i));
		%~ data1 = squeeze(data1(i,:,:));
		%~ 
		%~ data2 = RNX_data2.(const(const_i));
		%~ data2 = squeeze(data2(i,:,:));
		%~ 
		%~ if sum(sum(data1(2:end,:)));
		%~ 
		%~ off = 300*randn(size(data1,2),size(data_id,1));
		%~ 
		%~ 
		%~ for ep = 1:size(data1,2)
			%~ 
			%~ mjd = get_mjd_from_epoch(RNX_header1,ep);
			%~ ep2 = get_epoch_from_mjd(RNX_header2,mjd);
			%~ tic
			%~ off = 300*randn(size(data_id,1));
			%~ if(ep2~=0)
			%~ 
				%~ for var = 1:size(data_id,1)
					%~ 
					%~ data1(data_id(var,1)+1,ep) = data1(data_id(var,1)+1,ep) + off(ep,var);	
					%~ data2(data_id(var,2)+1,ep) = data2(data_id(var,2)+1,ep) + off(i);	
					%~ 
				%~ end
				%~ 
			%~ end
			%~ 
		%~ end
		%~ 
		%~ 
		%~ 
		%~ RNX_data1.(const(const_i))(i,:,:) = data1;
		%~ RNX_data2.(const(const_i))(i,:,:) = data2;
		%~ 
		%~ 
	%~ end
	%~ 
	%~ 
	%~ 
%~ end
%~ 
%~ toc
%~ end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
