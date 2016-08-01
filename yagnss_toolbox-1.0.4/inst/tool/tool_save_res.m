function [result_out] = tool_save_res(sta_pos,result,stat,process_time)
%%function [result_out] = tool_save_res(sta_pos,result,stat,process_time)
%%
%% Fill a structure array with results
%%
%% Clement Fontaine 2013-11-18
%%
%% Input : 
%% - X0 : Approximated values vector (initial)
%% - result : structure in output of calc_LS()
%% - stat : structure in output of calc_stat_indic()
%% - process_time : processing time
%%
%% Output : 
%% - result_out : structure array
%%		fields : 
%%          - sta_pos
%%			- X
%%			- Y
%%			- Z
%%			- E 
%%			- N
%%			- U
%%			- cdtr 
%%			- cGGTO 
%%			- cGPGL 
%%			- sigma02
%%			- n_iter
%%			- nb_GPS
%%			- nb_GLO
%%			- nb_GAL
%%			- V
%%			- Vnorm
%%			- index_sat
%%			- Qxx
%%			- conf_ell
%%			- Qenu
%%			- Corr_enu
%%			- GDOP
%%			- PDOP
%%			- HDOP
%%			- VDOP
%%			- TDOP
%%			- process_time
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if no computation
result_out.sta_pos = [0;0;0];
result_out.t = 0;
result_out.X = 0;
result_out.Y = 0;
result_out.Z = 0;

result_out.E = 0;
result_out.N = 0;
result_out.U = 0;

result_out.cdtr = 0;
result_out.cGGTO = 0;
result_out.cGPGL = 0;

result_out.sigma02 = 0;
result_out.n_iter = 0;

	
% sat number
result_out.nb_GPS = 0;
result_out.nb_GLO = 0;
result_out.nb_GAL = 0;
	
% Matrix
result_out.V = [];
result_out.Vnorm = [];
result_out.index_sat = cell(0);
result_out.Qxx = [];
	
% Quality factors
result_out.conf_ell = struct;
result_out.Qenu = [];
result_out.Corr_enu = [];
result_out.GDOP = 0;
result_out.PDOP = 0;
result_out.HDOP = 0;
result_out.VDOP = 0; 
result_out.TDOP =0;
	
result_out.process_time = 0;

if result.X ~=0


	% 2 cases : length(interm_results)==1 and length(interm_results)>1
	result_out.sta_pos = sta_pos(1:3);
	result_out.t = result.t;
	result_out.X = result.X;
	result_out.Y = result.Y;
	result_out.Z = result.Z;
	
	if(result.X~=0)
		[E,N,U] =  tool_cartloc_GRS80(sta_pos(1),sta_pos(2),sta_pos(3),result_out.X,result_out.Y,result_out.Z);
	else
		[E,N,U] = [0,0,0];	
	end
	
	result_out.E = E;
	result_out.N = N;
	result_out.U = U;
	
	% spp, DGNSS
	if(isfield(result,'cdtr'))
	
		result_out.cdtr = result.cdtr;
		result_out.cGGTO = result.cGGTO;
		result_out.cGPGL = result.cGPGL;
	
	% phase	
	else	
	
		result_out.cdtr = 0;
		result_out.cGGTO = 0;
		result_out.cGPGL = 0;
		
	end

	result_out.sigma02 = result.sigma02;
	result_out.n_iter = result.n_iter;
	
		
	% sat number
	result_out.nb_GPS = result.nb_GPS;
	result_out.nb_GLO = result.nb_GLO;
	result_out.nb_GAL = result.nb_GAL;
		
	% Matrix
	result_out.V = result.V;
	result_out.Vnorm = result.Vnorm;
	result_out.index_sat = result.sat_index;
	result_out.Qxx = result.Qxx;
		
	% Quality factors
	result_out.conf_ell = stat.ell;
	result_out.Qenu = stat.Qenu;
	result_out.Corr_enu = stat.Corr_enu;
	result_out.GDOP = stat.GDOP;
	result_out.PDOP = stat.PDOP;
	result_out.HDOP = stat.HDOP;
	result_out.VDOP = stat.VDOP;
	result_out.TDOP = stat.TDOP;
		
	result_out.process_time = process_time;
	
end
	
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
