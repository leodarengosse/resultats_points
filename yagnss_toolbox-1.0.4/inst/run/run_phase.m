function [result,interm_results,result_LS_code_base,result_LS_code_rover] = run_phase(rinex_o_base,rinex_o_rover,rinex_n,X0_base,options)
%% function [result,interm_results,result_LS_code_base,result_LS_code_rover] = run_phase(rinex_o_base,rinex_o_rover,rinex_n,X0_base,options)
%%
%% Estimate positions, dtr, time offset between constellations and ambiguities using double differences
%% Phase version : Double differences
%% Based on P3 and L3 obs
%%
%% Time offset : 
%% GGTO = GPS to Galileo Time Offset
%% GPGL = GPS to GLonass time offset
%% 
%% Clement Fontaine 2013-12-18
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
%%   nav : type of orbits
%%		   - 'brdc' : broadcasted ephemeris
%%         - 'sp3' : precise orbits
%%         default : 'brdc'    
%%   cut_off : elevation cut off in degree
%%         default : 3 degrees
%%   verbose : 1 = print information
%%		   default TODO: 1
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
%% - result : results of phase computation
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
%% - interm_results : intermediate results
%%
%%   {
%%       G_base
%%       G_rover
%%       calc
%%   }
%%
%% - result_LS_code_base : spp base computation results (same structure as result)
%% - result_LS_code_rover : spp rover computation results (same structure as result)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Close all figures
close('all')

%%%%% Default parameters
X0_rover = [0; 0; 0]; %[X, Y, Z] 
nav = 'brdc';
const = 'G';
sat_num = 32;
verbose=1;
cut_off = 3*pi/180;
Epoch_min = 1; % all file
N_epoch = -1; % all file
dir_out = '';

%%%%% Constant
c = 299792458.0;

%%%%% Options
if nargin==5
	[options,X0_rover,cut_off,verbose,nav,Epoch_min,N_epoch,const,sat_num,dir_out] = set_options(options,X0_rover,cut_off,verbose,nav,Epoch_min,N_epoch,const,sat_num,dir_out);
end

%%%%% Creation of output directory
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
tool_print_info('FUNCTION RUN_PHASE : DOUBLE DIFFERENCES ',1);
tool_print_info('----------------------------------------',1);

tool_print_info('',1);

tool_print_info(sprintf('BASE STATION : %s',base_name),1)
tool_print_info(sprintf('ROVER STATION : %s',rover_name),1)

tool_print_info('',1);

tool_print_info('OPTIONS : ',1);
tool_print_info(sprintf('\tCONSTELLATIONS : %s',const),1);
tool_print_info(sprintf('\tORBITS : %s',nav),1);
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

% Iono-free
eph_atx.const = 'G';
dE_base = 0;
dN_base = 0;
dU_base = 2.5*ATX_base_1.UP-1.5*ATX_base_2.UP;

% antenna height
dE_base = dE_base + RNX_header_base.dE;
dN_base = dN_base + RNX_header_base.dN;
dU_base = dU_base + RNX_header_base.dH;

[X0_base(1),X0_base(2),X0_base(3)] = corr_pos_atx(X0_base(1),X0_base(2),X0_base(3), dE_base, dN_base, dU_base, 1);

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

% Iono-free
eph_atx.const = 'G';
dE_rover = 0;
dN_rover = 0;
dU_rover = 2.5*ATX_rover_1.UP-1.5*ATX_rover_2.UP;

% antenna height
dE_rover = dE_rover + RNX_header_rover.dE;
dN_rover = dN_rover + RNX_header_rover.dN;
dU_rover = dU_rover + RNX_header_rover.dH;

[X0_rover(1),X0_rover(2),X0_rover(3)] = corr_pos_atx(X0_rover(1),X0_rover(2),X0_rover(3), dE_rover, dN_rover, dU_rover, 1);
                  
X0_rover = [X0_rover(1:3);0];                                        % initial cdtr

if (strcmp(nav,'brdc'))                              % initial cGGTO and cGPGL

	if isfield(NAV_header,'GPGA')
		cGGTO = c * NAV_header.GPGA(1);  % Approx GGTO, not the real formula which need mjd (cf official doc Galileo)
		X0_rover = [X0_rover;cGGTO];  
	else
		X0_rover = [X0_rover;0];  
	end

	if isfield(NAV_header,'GLGP')
		cGPGL = c * NAV_header.GLGP(1); % in rinex, GPGL is not defined, but  - GLGP is present
		X0_rover = [X0_rover;cGPGL];  
	else
		X0_rover = [X0_rover;0];  
	end
	
else % sp3
	X0_rover = [X0_rover;0;0]; 
end

%%%%% Set options

comp_options.X0 = X0_rover;
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
result = struct;
interm_results = struct;
result_LS_code_base=struct;
result_LS_code_rover=struct;

global_process_time = now;

			
%%%%% first step : Preprocessing

tool_print_info('PREPROCESSING',1)
tool_print_info('',1);

%%%%% matrix initialisation
t = zeros(sat_num*N_epoch,1);
PosSat = zeros(sat_num*N_epoch,6);
Dobs = zeros(sat_num*N_epoch,2);
PobsL1 = zeros(sat_num*N_epoch,2);
PobsL2 = zeros(sat_num*N_epoch,2);
PobsL3 = zeros(sat_num*N_epoch,2);
AzSat = zeros(sat_num*N_epoch,2);
ElevSat = zeros(sat_num*N_epoch,2);
Dtropo = zeros(sat_num*N_epoch,2);
sat_index = cell(sat_num*N_epoch,2);



N_val = 0;

mean_t = 0;

N_epoch = comp_options.N_epoch;
Epoch_min = comp_options.Epoch_min;

for i = 1:N_epoch

	tool_print_info('',0)
	tool_print_info(sprintf('EPOCH %d/%d : ',i,N_epoch),1)
	tool_print_info(sprintf('\tEST REM TIME : %2d MIN %0d S   ',floor((N_epoch-i)*mean_t/60),fix(((N_epoch-i)*mean_t/60-floor((N_epoch-i)*mean_t/60))*60)),0);
	
	process_time = now;

	% preprocessing for epoch i
	[G1,G2,calc,result_LS_code_base_i,result_LS_code_rover_i]=calc_preprocessing_phase(RNX_header_base, RNX_data_base, RNX_header_rover, RNX_data_rover, NAV_header, NAV_data, X0_base, Epoch_min+i-1, comp_options);
	
	tool_print_info(sprintf('\tSAT : %2d GPS, %2d GLO, %2d GAL',calc.nb_GPS,calc.nb_GLO,calc.nb_GAL),0);

	process_time = (now - process_time)*86400;

	% save data for phase positionning
	nb_loc_val = size(calc.PosSat,1);
	
	t(N_val+1:N_val + nb_loc_val,:) = calc.t.mjd*ones(nb_loc_val,1);

	PosSat(N_val+1:N_val + nb_loc_val,:) = calc.PosSat;
	Dobs(N_val+1:N_val + nb_loc_val,:) = calc.Dobs;
	PobsL1(N_val+1:N_val + nb_loc_val,:) = calc.PobsL1;
	PobsL2(N_val+1:N_val + nb_loc_val,:) = calc.PobsL2;
	PobsL3(N_val+1:N_val + nb_loc_val,:) = calc.PobsL3;
	AzSat(N_val+1:N_val + nb_loc_val,:) = calc.AzSat;
	ElevSat(N_val+1:N_val + nb_loc_val,:) = calc.ElevSat;
	Dtropo(N_val+1:N_val + nb_loc_val,:) = calc.Dtropo;
	
	sat_index(N_val+1:N_val + nb_loc_val,1) = calc.sat_index(:,1);
	sat_index(N_val+1:N_val + nb_loc_val,2) = calc.sat_index(:,2);

	N_val = N_val + nb_loc_val;
	
	interm_results(i).G_base = G1;
	interm_results(i).G_rover = G2;
	interm_results(i).calc = calc;
	
	% save LS code
	
	% Statistic indicators computation on code
	% base
	[stat_code_base]=calc_stat_indic(result_LS_code_base_i.Qxx,result_LS_code_base_i.X,result_LS_code_base_i.Y,result_LS_code_base_i.Z);
	% rover
	[stat_code_rover]=calc_stat_indic(result_LS_code_rover_i.Qxx,result_LS_code_rover_i.X,result_LS_code_rover_i.Y,result_LS_code_rover_i.Z);
	
	% phase center and antenna height correction
	% base
	[result_LS_code_base_i.X,result_LS_code_base_i.Y,result_LS_code_base_i.Z] = corr_pos_atx(result_LS_code_base_i.X,result_LS_code_base_i.Y,result_LS_code_base_i.Z, dE_base, dN_base, dU_base, 2);

	% rover
	[result_LS_code_rover_i.X,result_LS_code_rover_i.Y,result_LS_code_rover_i.Z] = corr_pos_atx(result_LS_code_rover_i.X,result_LS_code_rover_i.Y,result_LS_code_rover_i.Z, dE_rover, dN_rover, dU_rover,2);
	
	% save results and -> for matlab, need to preallocate result 
	% otherwise error 'the input character is not valid in MATLAB statements or expressions'
	
	% base
	if i==1
		result_LS_code_base = tool_save_res(base_sta_pos,result_LS_code_base_i,stat_code_base,process_time);	
	else
		result_LS_code_base(i) = tool_save_res(base_sta_pos,result_LS_code_base_i,stat_code_base,process_time);
	end
	
	% rover
	if i==1
		result_LS_code_rover = tool_save_res(rover_sta_pos,result_LS_code_rover_i,stat_code_rover,process_time);	
	else
		result_LS_code_rover(i) = tool_save_res(rover_sta_pos,result_LS_code_rover_i,stat_code_rover,process_time);
	end
	
	mean_t = mean([process_time mean_t]);

end

%%%%% redim matrix
t = t(1:N_val,:);
PosSat = PosSat(1:N_val,:);
Dobs = Dobs(1:N_val,:);
PobsL1 = PobsL1(1:N_val,:);
PobsL2 = PobsL2(1:N_val,:);
PobsL3 = PobsL3(1:N_val,:);
AzSat = AzSat(1:N_val,:);
ElevSat = ElevSat(1:N_val,:);
Dtropo = Dtropo(1:N_val,:);
sat_index = sat_index(1:N_val,:);

%%%%% Cycle slip 

tool_print_info('',1);
tool_print_info('CYCLE SLIP DETECTION',1)
tool_print_info('',1);

[t,Pobs,Dobs,PosSat,ElevSat,AzSat,Dtropo,sat_index,amb_index,amb0] = calc_cycle_slip(t,PobsL1,PobsL2,PobsL3,Dobs,PosSat,ElevSat,AzSat,Dtropo,sat_index);

%%%%% Select pivot satellite

tool_print_info('',1);
tool_print_info('PIVOT SATELLITES SELECTION',1)
tool_print_info('',1);

[pivot,interv] = calc_pivot(t,ElevSat,sat_index);

%%%%% Least Squares

tool_print_info('',1);
tool_print_info('LS ESTIMATION',1)
tool_print_info('',1);

[result_LS] = calc_LS_phase(t,Pobs,PosSat,ElevSat,Dtropo,sat_index,amb_index,amb0,X0_base,X0_rover,pivot,interv);

%%%%% Phase center and antenna height correction
[result_LS.X,result_LS.Y,result_LS.Z] = corr_pos_atx(result_LS.X,result_LS.Y,result_LS.Z, dE_rover, dN_rover, dU_rover,2);

%%%%% Statistic indicators computation
[stat]=calc_stat_indic(result_LS.Qxx,result_LS.X,result_LS.Y,result_LS.Z);

global_process_time = (now - global_process_time)*86400;


%%%%% save results
[result] = tool_save_res(rover_sta_pos,result_LS,stat,global_process_time);


%%%%% plot - spp

if ~strcmp(dir_out,'')
	save_plot_result_base = strcat(dir_out,filesep(),'base_spp');
	save_plot_result_rover = strcat(dir_out,filesep(),'rover_spp');
else
		save_plot_result_base = '';
		save_plot_result_rover = '';
end	

% plot results base
plot_results_code(result_LS_code_base,rinex_o_base(end-11:end-8),save_plot_result_base);
% plot results rover
plot_results_code(result_LS_code_rover,rinex_o_rover(end-11:end-8),save_plot_result_rover);


%%%%% plot - phase
% plot histo - residus
if ~strcmp(dir_out,'')
	save_hist_res_norm = strcat(dir_out,filesep(),'hist_res_norm');
else
	save_hist_res_norm = '';
end
figure()
plot_hist(result.Vnorm,40,'Histogram : normalized residuals','Value',save_hist_res_norm);

%~ if ~strcmp(dir_out,'')
	%~ save_hist_res = strcat(dir_out,filesep(),'hist_res');
%~ else
	%~ save_hist_res = '';
%~ end
%~ figure()
%~ plot_hist(result.V,40,'Histogram : residuals','m',save_hist_res);

% correlation matrix
if ~strcmp(dir_out,'')
	save_corr = strcat(dir_out,filesep(),'corr');
else
	save_corr = '';
end	
figure()
plot_corr_matrix(result.Corr_enu,save_corr);
	
% skyplot
if ~strcmp(dir_out,'')
	save_skyplot = strcat(dir_out,filesep(),'skyplot_base');
else
	save_skyplot = '';
end	
mat_skyplot = [t,AzSat(:,1),ElevSat(:,1)];
plot_skyplot(mat_skyplot,sat_index,X0_base(1),X0_base(2),X0_base(3),rinex_o_base(end-11:end-8),save_skyplot);

if ~strcmp(dir_out,'')
	save_skyplot = strcat(dir_out,filesep(),'skyplot_rover');
else
	save_skyplot = '';
end	
mat_skyplot = [t,AzSat(:,2),ElevSat(:,2)];
plot_skyplot(mat_skyplot,sat_index,result.X,result.Y,result.Z,rinex_o_rover(end-11:end-8),save_skyplot);

tool_print_info('',1);

%%%%% report
date_now = datestr(now,'yy_mm_dd_HHMMss');

if ~strcmp(dir_out,'')
	save_report_phase = strcat(dir_out,filesep(),'results_phase_',date_now,'.txt');
	save_report_code_base = strcat(dir_out,filesep(),'results_code_base_',date_now,'.txt');
	save_report_code_rover = strcat(dir_out,filesep(),'results_code_rover_',date_now,'.txt');
else
	save_report_phase = strcat('results_phase',date_now,'.txt');
	save_report_code_base = strcat('results_code_base',date_now,'.txt');
	save_report_code_rover = strcat('results_code_rover',date_now,'.txt');
end

tool_print_info('',1);
tool_print_info('----------------------------------------',1);
tool_print_info('CODE : BASE',1);	
tool_print_info('----------------------------------------',1);
tool_print_info('',1);

tool_save_report(result_LS_code_base,save_report_code_base);

tool_print_info('',1);
tool_print_info('----------------------------------------',1);
tool_print_info('CODE : ROVER',1);	
tool_print_info('----------------------------------------',1);
tool_print_info('',1);
	
tool_save_report(result_LS_code_rover,save_report_code_rover);	

tool_print_info('',1);
tool_print_info('----------------------------------------',1);
tool_print_info('PHASE',1);	
tool_print_info('----------------------------------------',1);
tool_print_info('',1);

tool_save_report(result,save_report_phase);	


tool_print_info('',1);
tool_print_info(sprintf('TOTAL PROCESSING TIME : %0.3f S',global_process_time),1);	
tool_print_info('',1);

tool_print_info('----------------------------------------',1);
tool_print_info('         END FUNCTION RUN_PHASE',1);
tool_print_info('----------------------------------------',1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [options,X0,cut_off,verbose,nav,Epoch_min,N_epoch,const,sat_num,dir_out] = set_options(options,X0,cut_off,verbose,nav,Epoch_min,N_epoch,const,sat_num,dir_out)
%% function [options,X0,cut_off,verbose,nav,Epoch_min,N_epoch,const,sat_num,dir_out] = set_options(options,X0,cut_off,verbose,nav,Epoch_min,N_epoch,const,sat_num,dir_out)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if isfield(options,'X0')
        X0 = options.X0(:); % in column
    else
		options.X0 = X0;
    end  
           
    if isfield(options,'cut_off')
        cut_off = options.cut_off*pi/180; % input in degree, transformation to radian
    else    
		options.cut_off = cut_off;
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
