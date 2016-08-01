function [result_base,interm_results_base,result_rover,interm_results_rover] = run_DGNSS(rinex_o_base,rinex_o_rover,rinex_n,X0_base,options)
%% function [result_base,interm_results_base,result_rover,interm_results_rover] = run_DGNSS(rinex_o_base,rinex_o_rover,rinex_n,X0_base,options)
%%
%% Estimate positions, dtr and time offset between constellations using 
%% Pseudo-range. 
%% DGNSS version : Differential GNSS
%%
%% Time offset : 
%% GGTO = GPS to Galileo Time Offset
%% GPGL = GPS to GLonass time offset
%% 
%% Clement Fontaine 2013-11-14
%%
%% Input : 
%% - rinex_o_base : observation RINEX name of base station
%% - rinex_o_rover : observation RINEX name of rover station
%% - rinex_n : navigation RINEX name OR sp3 name cell array {GPS.sp3;GLO.sp3;GAL.sp3}
%% - X0_base : base station coordinates [X;Y;Z] in m
%% - options : structure containing options of computation :
%% 
%%   options = 
%% 
%% {
%%   X0 : approximated coordinates of rover station (column vector of 3 elements [X;Y;Z]) 
%%        default : X0 = [0 ;0; 0];
%%   const : constellations used ('G' = GPS, 'R' = Glonass, 'E' = Galileo, for multi-constellation concatenate chars)
%%        default : 'G'
%%   freq : type of used data [C1,P2,iono_free] 
%%        default : iono_free
%%         - 'F1' : use F1 frequency obs
%%         - 'F2' : use F2 frequency obs (or F5 for Galileo)
%%         - 'iono_free' : ionosphere-free combination
%%		  default : 'iono_free'
%%   iono : type of correction :
%%		   - 'klobuchar' : klobuchar modelization -> if  nav = 'brdc'
%%		   - 'none' : no correction (if freq = 'iono_free', iono set to 'none')
%%         default : 'none' 
%%   nav : type of orbits
%%		   - 'brdc' : broadcasted ephemeris
%%         - 'sp3' : precise orbits
%%         default : 'brdc'    
%%   cut_off : elevation cut off in degree
%%         default : 3 degrees
%%   verbose : 1 = print information
%%		   default : 1
%%   Epoch_min : number of first RINEX observation to use
%%         default : 1
%%   N_epoch : number of epochs to use. if < 0 or > size(RNX_data.G,3), 
%%         set to size(RNX_data.G,3)
%%         default : -1
%%   dir_out : output directory
%%		   default : ''
%% }
%%
%% Output : 
%% - result_base : result structure array for base station
%%
%%   {   
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
%%   }
%%
%% - interm_results_base : intermediate results for base station
%%
%%   {
%%     1x20 struct array containing the fields:
%%
%%       G
%%       calc
%%   }
%% - result_rover : results for rover station (same structure than result_base)
%% - interm_results_rover : intermediate results for rover station (same structure than interm_results_base)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Close all figures
close('all')

%%%%% Default parameters
X0_rover = [0; 0; 0]; %[X, Y, Z] 
freq = 'iono_free';
iono = 'none'; % iono_free -> no more correction
nav = 'brdc';
const = 'G';
sat_num = 32;
verbose=1;
cut_off = 3*pi/180;
Epoch_min = 1; % Beginning
N_epoch = -1; % all file
dir_out = '';

%%%%% Constant
c = 299792458.0;

%%%%% Options
if nargin==5

	[options,X0_rover,freq,cut_off,iono,verbose,nav,Epoch_min,N_epoch,const,sat_num,dir_out] = set_options(options,X0_rover,freq,cut_off,iono,verbose,nav,Epoch_min,N_epoch,const,sat_num,dir_out);

end

%%%%% Create output directory
tool_create_dirs(dir_out);

%%%%% Log file and terminal output initialization
if strcmp(dir_out,'')
	log_name = strcat('.',filesep(),'log_',datestr(now,'yy_mm_dd_HHMMss'),'.txt');
else
	log_name = strcat(dir_out,filesep(),'log_',datestr(now,'yy_mm_dd_HHMMss'),'.txt');
end

global LOG_FILE; 
LOG_FILE = log_name;

global VERBOSE;
VERBOSE = verbose;

%%%%% Displaying parameters
if length(rinex_o_base)<12
	rinex_o_base = strcat(rinex_o_base,repmat('x',1,12-length(rinex_o)))
end

if length(rinex_o_rover)<12
	rinex_o_rover = strcat(rinex_o_rover,repmat('x',1,12-length(rinex_o)))
end

base_name = rinex_o_base(end-11:end-8);
rover_name = rinex_o_rover(end-11:end-8);

tool_print_info('----------------------------------------',1)
tool_print_info('FUNCTION RUN_DGNSS : DIFFERENTIAL GNSS',1);
tool_print_info('----------------------------------------',1);

tool_print_info('',1);

tool_print_info(sprintf('BASE STATION : %s',base_name),1)
tool_print_info(sprintf('ROVER STATION : %s',rover_name),1)

tool_print_info('',1);

tool_print_info('OPTIONS : ',1);
tool_print_info(sprintf('\tCONSTELLATIONS : %s',const),1);
tool_print_info(sprintf('\tORBITS : %s',nav),1);
tool_print_info(sprintf('\tFREQUENCY : %s',freq),1);
tool_print_info(sprintf('\tIONO : %s',iono),1);
tool_print_info(sprintf('\tCUT OFF (DEG): %d',cut_off*180/pi),1);


tool_print_info('',1);
tool_print_info('DATA LOADING',1);
tool_print_info('',1);
	

%%%%% Data loading
% --> Base
[RNX_header_base,RNX_data_base]=load_rinex_o(rinex_o_base);
if RNX_header_base.VERSION==0
	return;
end
% --> Rover
[RNX_header_rover,RNX_data_rover]=load_rinex_o(rinex_o_rover);
if RNX_header_rover.VERSION==0
	return;
end
% --> Nav or sp3
if strcmp(nav,'brdc')                                  % NAV
	[NAV_header,NAV_data]=load_rinex_n(rinex_n);
	if NAV_header.VERSION==0
		return;
	end
else                                                   % SP3
	[NAV_header,NAV_data]=load_sp3(rinex_n);
	if (~isfield(NAV_header.G,'Version') && ~isfield(NAV_header.R,'Version') && ~isfield(NAV_header.E,'Version'))
		return;
	end
end

% Add SA
%~ [RNX_header_base,RNX_data_base,RNX_header_rover,RNX_data_rover] = tool_add_SA(RNX_header_base,RNX_data_base,RNX_header_rover,RNX_data_rover);

%%%%% Set epoch min and epoch max
% Epoch -> select data range

% For rover
if(Epoch_min < 1 || Epoch_min > size(RNX_data_rover.G,3))
	Epoch_min_rover = 1;
else
	Epoch_min_rover = Epoch_min;
end

if(N_epoch < 1 || N_epoch > size(RNX_data_rover.G,3)-Epoch_min_rover+1)
	N_epoch_rover = size(RNX_data_rover.G,3)-Epoch_min_rover+1;
else
	N_epoch_rover = N_epoch;
end

% tmin and tmax for rover
tmin_rover = get_mjd_from_epoch(RNX_header_rover,Epoch_min_rover);
tmax_rover = get_mjd_from_epoch(RNX_header_rover,Epoch_min_rover + N_epoch_rover - 1);

% For base
Epoch_min_base = get_epoch_from_mjd(RNX_header_base,tmin_rover,86400);
Epoch_max_base = get_epoch_from_mjd(RNX_header_base,tmax_rover,86400);
N_epoch_base = Epoch_max_base - Epoch_min_base + 1;

% test corresponding data
if (tmin_rover == 0 || tmax_rover == 0 || Epoch_min_base == 0 || Epoch_max_base == 0)

	tool_print_info('ERROR : No corresponding obs between Base and Rover)',3);
	return;	

end

%%%%% Initial values -> base
X0_base = [X0_base(1:3);0;0;0];

% base station position 
base_sta_pos = X0_base(1:3);

% phase center and antenna height correction

% calibration -> base
antex_file = 'igs08.atx';
[ATX_header,ATX_data] = load_antex(antex_file);
[ATX_base_1] = get_antex(ATX_header, ATX_data, RNX_header_base.ANT_TYPE, 'G01');
[ATX_base_2] = get_antex(ATX_header, ATX_data, RNX_header_base.ANT_TYPE, 'G02');

if strcmp(freq, 'F1')            % First frequency
	dE_base = 0;
	dN_base = 0;
	dU_base = ATX_base_1.UP;
elseif strcmp(freq, 'F2')        % Second frequency
	dE_base = 0;
	dN_base = 0;
	dU_base = ATX_base_2.UP;
else                             % Iono-free
	eph_atx.const = 'G';
	dE_base = 0;
	dN_base = 0;
	dU_base = 2.5*ATX_base_1.UP-1.5*ATX_base_2.UP;
end

% antenna height
dE_base = dE_base + RNX_header_base.dE;
dN_base = dN_base + RNX_header_base.dN;
dU_base = dU_base + RNX_header_base.dH;

% station position to phase center (approximate coordinates)
[X0_base(1),X0_base(2),X0_base(3)] = corr_pos_atx(X0_base(1),X0_base(2),X0_base(3), dE_base, dN_base, dU_base,1);

%%%%% Initial values -> rover   
if(sum(X0_rover)==0)
	X0_rover = [RNX_header_rover.X;RNX_header_rover.Y;RNX_header_rover.Z]; % initial position
end      

% rover station position 
rover_sta_pos = X0_rover(1:3);

% phase center and antenna height correction
% calibration -> rover
[ATX_rover_1] = get_antex(ATX_header, ATX_data, RNX_header_rover.ANT_TYPE, 'G01');
[ATX_rover_2] = get_antex(ATX_header, ATX_data, RNX_header_rover.ANT_TYPE, 'G02');

if strcmp(freq, 'F1')            % First frequency
	dE_rover = 0;
	dN_rover = 0;
	dU_rover = ATX_rover_1.UP;
elseif strcmp(freq, 'F2')        % Second frequency
	dE_rover = 0;
	dN_rover = 0;
	dU_rover = ATX_rover_2.UP;
else                             % Iono-free
	eph_atx.const = 'G';
	dE_rover = 0;
	dN_rover = 0;
	dU_rover = 2.5*ATX_rover_1.UP-1.5*ATX_rover_2.UP;
end

% antenna height
dE_rover = dE_rover + RNX_header_rover.dE;
dN_rover = dN_rover + RNX_header_rover.dN;
dU_rover = dU_rover + RNX_header_rover.dH;

[X0_rover(1),X0_rover(2),X0_rover(3)] = corr_pos_atx(X0_rover(1),X0_rover(2),X0_rover(3), dE_rover, dN_rover, dU_rover,1);
                        
X0_rover = [X0_rover(1:3);0];                                        % initial cdtr
X0_base = [X0_base(1:3);0];  

if (strcmp(nav,'brdc'))                              % initial cGGTO and cGPGL

	if isfield(NAV_header,'GPGA')
		cGGTO = c * NAV_header.GPGA(1);  % Approx GGTO, not the real formula which need mjd (cf official doc Galileo)
		X0_rover = [X0_rover;cGGTO];  
		X0_base = [X0_base;cGGTO];  
	else
		X0_rover = [X0_rover;0];  
		X0_base = [X0_base;0];  
	end

	if isfield(NAV_header,'GLGP')
		cGPGL = c * NAV_header.GLGP(1); % in rinex, GPGL is not defined, but  - GLGP is present
		X0_rover = [X0_rover;cGPGL];  
		X0_base = [X0_base;cGPGL];  
	else
		X0_rover = [X0_rover;0];  
		X0_base = [X0_base;0];  
	end
	
else % sp3
	X0_rover = [X0_rover;0;0]; 
	X0_base = [X0_base;0;0]; 
end

%%%%% Set options 

comp_options.X0 = X0_rover;
comp_options.freq = freq;
comp_options.iono = iono;
comp_options.nav = nav;
comp_options.const = const;
comp_options.sat_num = sat_num;
comp_options.verbose = verbose;
comp_options.cut_off = cut_off;
comp_options.Epoch_min = Epoch_min_rover;
comp_options.N_epoch = N_epoch_rover;


tool_print_info('',1);
tool_print_info(sprintf('BASE COORDINATES (M): %0.3f %0.3f %0.3f',X0_base(1:3)),1);
tool_print_info('',1);
tool_print_info('INITIAL VALUES',1);
tool_print_info(sprintf('\tINITIAL ROVER COORDINATES (M): %0.3f %0.3f %0.3f',X0_rover(1:3)),1);
tool_print_info(sprintf('\tINITIAL DTR (S) : %0.3f', X0_rover(4)),1);
tool_print_info(sprintf('\tINITIAL cGGTO (M) : %0.3f', X0_rover(5)),1);
tool_print_info(sprintf('\tINITIAL cGPGL (M) : %0.3f', X0_rover(6)),1);
tool_print_info('',1);
	
%%%%% Output
result_base = struct;
interm_results_rover = struct;
result_base = struct;
interm_results_rover = struct;


global_process_time = now;

		
mean_t = 0;

N_epoch = comp_options.N_epoch;
Epoch_min = comp_options.Epoch_min;

for i = 1:N_epoch
	
	tool_print_info('',0)
	tool_print_info(sprintf('EPOCH %d/%d : ',i,N_epoch),1)
	tool_print_info(sprintf('\tEST REM TIME : %d MIN %d S   ',floor((N_epoch-i)*mean_t/60),fix(((N_epoch-i)*mean_t/60-floor((N_epoch-i)*mean_t/60))*60)),0);
	
	process_time = now;

	% Preprocessing
	[G_base, G_rover, calc, result_LS_base]=calc_preprocessing_DGNSS(RNX_header_base, RNX_data_base, RNX_header_rover, RNX_data_rover, NAV_header, NAV_data, X0_base, Epoch_min + i - 1, comp_options);

	tool_print_info(sprintf('\tSAT : %2d GPS, %2d GLO, %2d GAL',calc.nb_GPS,calc.nb_GLO,calc.nb_GAL),0);
				
	t = calc.t.mjd*ones(calc.nb_sat,1);
	PosSat = calc.PosSat(:,4:6);
	Dobs = calc.Dobs(:,2);
	ElevSat = calc.ElevSat(:,2);
	sat_index = calc.sat_index;

	% Least Square computation
	[result_LS_rover] = calc_LS_code(t,PosSat,Dobs,ElevSat,sat_index,X0_rover,comp_options);
	
	% Statistic indicators computation
	% base
	[stat_base]=calc_stat_indic(result_LS_base.Qxx,result_LS_base.X,result_LS_base.Y,result_LS_base.Z);
	% rover
	[stat_rover]=calc_stat_indic(result_LS_rover.Qxx,result_LS_rover.X,result_LS_rover.Y,result_LS_rover.Z);
	
	% processing time
	process_time = (now - process_time)*86400;
	mean_t = mean([process_time mean_t]);
	
	
	% update initial coordinates + cdtr-> calc_LS_code output = input if no estimation
	X0_rover(1:4) = [result_LS_rover.X;result_LS_rover.Y;result_LS_rover.Z;result_LS_rover.cdtr];
	options.X0 = X0_rover;
		
	
	% phase center and antenna height correction
	% base
	[result_LS_base.X, result_LS_base.Y, result_LS_base.Z] = corr_pos_atx(result_LS_base.X, result_LS_base.Y, result_LS_base.Z, dE_base, dN_base, dU_base,2);

	% rover
	[result_LS_rover.X, result_LS_rover.Y, result_LS_rover.Z] = corr_pos_atx(result_LS_rover.X, result_LS_rover.Y, result_LS_rover.Z, dE_rover, dN_rover, dU_rover,2);
	
	% save results and intermediates values -> for matlab, need to preallocate result 
	% otherwise error 'the inpur character is not valid in MATLAB statements or expressions
	
	% base
	if i==1
		result_base = tool_save_res(base_sta_pos,result_LS_base,stat_base,process_time);	
		interm_results_base.G = G_base;
		interm_results_base.calc = calc;
	else
		result_base(i) = tool_save_res(base_sta_pos,result_LS_base,stat_base,process_time);
		interm_results_base(i).G = G_base;
		interm_results_base(i).calc = calc;
	end
	
	% rover
	if i==1
		result_rover = tool_save_res(rover_sta_pos,result_LS_rover,stat_rover,process_time);	
		interm_results_rover.G = G_rover;
		interm_results_rover.calc = calc;
	else
		result_rover(i) = tool_save_res(rover_sta_pos,result_LS_rover,stat_rover,process_time);
		interm_results_rover(i).G = G_rover;
		interm_results_rover(i).calc = calc;
	end
	

			
end

tool_print_info('',1);

% plots

if ~strcmp(dir_out,'')
	save_plot_result_base = strcat(dir_out,filesep(),'base');
	save_plot_result_rover = strcat(dir_out,filesep(),'rover');
else
		save_plot_result_base = '';
		save_plot_result_rover = '';
end	

% plot results base
plot_results_code(result_base,rinex_o_base(end-11:end-8),save_plot_result_base);
% plot results rover
plot_results_code(result_rover,rinex_o_rover(end-11:end-8),save_plot_result_rover);


% skyplot

if ~strcmp(dir_out,'')
	save_skyplot_base = strcat(dir_out,filesep(),'base',filesep(),'skyplot');
else
		save_skyplot_base = '';
end	

[mat_skyplot_base,sat_index_base] = tool_prep_skyplot(interm_results_base);	
plot_skyplot(mat_skyplot_base,sat_index_base,base_sta_pos(1),base_sta_pos(2),base_sta_pos(3),rinex_o_base(end-11:end-8),save_skyplot_base);


global_process_time = (now - global_process_time)*86400;

tool_print_info('',1);
tool_print_info('BASE',1);
tool_print_info('',1);

% report_base
if ~strcmp(dir_out,'')
	save_report = strcat(dir_out,filesep(),'results_base_',datestr(now,'yy_mm_dd_HHMMss'),'.txt');
else
	save_report = strcat('results_base_',datestr(now,'yy_mm_dd_HHMMss'),'.txt');
end
tool_save_report(result_base,save_report);


tool_print_info('',1);
tool_print_info('ROVER',1);
tool_print_info('',1);

% report rover
if ~strcmp(dir_out,'')
	save_report = strcat(dir_out,filesep(),'results_rover_',datestr(now,'yy_mm_dd_HHMMss'),'.txt');
else
	save_report = strcat('results_rover_',datestr(now,'yy_mm_dd_HHMMss'),'.txt');
end
tool_save_report(result_rover,save_report);


tool_print_info('',1);
tool_print_info(sprintf('TOTAL PROCESSING TIME : %0.3f S',global_process_time),1);	
tool_print_info('',1);

tool_print_info('----------------------------------------',1);
tool_print_info('         END FUNCTION RUN_DGNSS',1);
tool_print_info('----------------------------------------',1);
	
%~ fclose('all');

end





function [options,X0,freq,cut_off,iono,verbose,nav,Epoch_min,N_epoch,const,sat_num,dir_out] = set_options(options,X0,freq,cut_off,iono,verbose,nav,Epoch_min,N_epoch,const,sat_num,dir_out)
%% function [options,X0,freq,cut_off,iono,verbose,nav,global_pos,Epoch_min,N_epoch,const,sat_num,dir_out] = set_options(options,X0,freq,cut_off,iono,verbose,nav,global_pos,Epoch_min,N_epoch,const,sat_num,dir_out)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if isfield(options,'X0')
        X0 = options.X0(:); % in column
    else
		options.X0 = X0;
	end  
         
    if isfield(options,'freq')
        if (strcmp(options.freq,'F1') || strcmp(options.freq,'F2') || strcmp(options.freq,'iono_free')) 
            freq = options.freq;  
        else    
			options.freq = 'iono_free';
        end 
    else    
		options.freq = 'iono_free';
    end 
  
	if isfield(options,'cut_off')
        cut_off = options.cut_off*pi/180; % input in degree, transformation to radian
    else    
		options.cut_off = cut_off;
    end
    
    if isfield(options,'iono') 
        if (strcmp(options.iono,'klobuchar') || strcmp(options.iono,'none')) 
            iono = options.iono; 
        else 
			options.iono = iono;
        end
    else 
		options.iono = iono;
    end  
    
    if isfield(options,'verbose') 
        
        if options.verbose==0
            verbose = options.verbose;   
        else
			options.verbose = verbose;
        end 
    else
		options.verbose = verbose;
    end  
    
    if isfield(options,'nav') 
        if (strcmp(options.nav,'sp3'))
            nav = 'sp3'; 
        else 
			options.nav = nav;
        end  
    else 
		options.nav = nav;
    end  
    
    if isfield(options,'const')
		sat_num = 0;
		const = '';
		if regexp(options.const,'G')>0
			const = strcat(const,'G');
			sat_num = sat_num + 32;
		end
		if regexp(options.const,'R')>0
			const = strcat(const,'R');
			sat_num = sat_num + 32;
		end	
		if regexp(options.const,'E')>0
			const = strcat(const,'E');
			sat_num = sat_num + 32;
		end      
    else 
		options.const = const;
    end  
    
    if isfield(options,'Epoch_min') 
        if (options.Epoch_min>1)
            Epoch_min = options.Epoch_min;
        else 
			options.Epoch_min = Epoch_min;
        end  
    else 
		options.Epoch_min = Epoch_min;
    end 
     
    if isfield(options,'N_epoch') 
        if (options.N_epoch>0)
            N_epoch = options.N_epoch; 
        else 
			options.N_epoch = N_epoch;   
        end  
    else 
		options.N_epoch = N_epoch;
    end  
     
    if isfield(options,'dir_out') 
        dir_out = options.dir_out;
    else 
		options.dir_out = dir_out;
    end  

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
