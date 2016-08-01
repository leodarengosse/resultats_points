function [result,interm_results] = run_spp(rinex_o,rinex_n,options)
%% function [result,interm_results] = run_spp(rinex_o,rinex_n,options)
%%
%% Estimate positions, dtr and time offset between constellations using 
%% Pseudo-range. 
%% 
%% Time offset : 
%% GGTO = GPS to Galileo Time Offset
%% GPGL = GPS to GLonass time offset
%%
%% Clement Fontaine 2013-11-14
%%
%% Input : 
%% - rinex_o : observation RINEX name
%% - rinex_n : navigation RINEX name OR sp3 name cell array {GPS.sp3;GLO.sp3;GAL.sp3}
%% - options : structure containing options of computation :
%% 
%%   options = 
%% 
%% {
%%   X0 : approximated coordinates (column vector of 3 elements [X;Y;Z]) 
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
%%   global_pos : 1 for one pos for the whole RINEX, 0 for one pos per epoch
%%         default : 0
%%   dir_out : output directory
%%		   default : ''
%% }
%%
%% Output : 
%% - result : result structure array
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
%%     1x20 struct array containing the fields:
%%
%%       G
%%       calc
%%   }
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = 299792458.0;

%%%%% Close all figures
close('all')

%%%%% Default parameters
X0 = [0; 0; 0]; %[X, Y, Z] 
freq = 'iono_free';
iono = 'none'; % iono_free -> no more correction
nav = 'brdc';
const = 'G';
sat_num = 32;
verbose=1;
cut_off = 3*pi/180;
global_pos = 0; % one position, cGGTO and cGPGL per RNX, 1 cdtr per epoch
Epoch_min = 1; % all file
N_epoch = -1; % all file
dir_out = '';

%%%%% Output preallocation
result = struct;
interm_results = struct; 

%%%%% Options
% Get options if they are specified
if nargin==3

	[options,X0,freq,cut_off,iono,verbose,nav,global_pos,Epoch_min,N_epoch,const,sat_num,dir_out] = set_options(options,X0,freq,cut_off,iono,verbose,nav,global_pos,Epoch_min,N_epoch,const,sat_num,dir_out);
      
end

%%%%% Create output directory
tool_create_dirs(dir_out);

%%%%% Initialize log file and terminal output
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
if length(rinex_o)<12
	rinex_o = strcat(rinex_o,repmat('x',1,12-length(rinex_o)))
end

rinex_name = rinex_o(end-11:end-8);

tool_print_info('----------------------------------------',1)
tool_print_info('FUNCTION RUN_SPP : SINGLE POINT POSITION',1);
tool_print_info('----------------------------------------',1);

tool_print_info('',1);

tool_print_info(sprintf('STATION : %s',rinex_name),1)

tool_print_info('',1);

tool_print_info('OPTIONS : ',1);
tool_print_info(sprintf('\tCONSTELLATIONS : %s',const),1);
tool_print_info(sprintf('\tORBITS : %s',nav),1);
tool_print_info(sprintf('\tFREQUENCY : %s',freq),1);
tool_print_info(sprintf('\tIONO : %s',iono),1);
tool_print_info(sprintf('\tCUT OFF (DEG): %d',cut_off*180/pi),1);

if global_pos==1
	tool_print_info(sprintf('\tCompute one POS, GGTO and GPGL for all rinex, and 1 dtr per epoch'),1);
else
	tool_print_info(sprintf('\tCompute one POS, dtr, GGTO and GPGL per epoch'),1);
end

tool_print_info('',1);
tool_print_info('DATA LOADING',1);
tool_print_info('',1);
		

%%%%% Data loading
% --> Obs
[RNX_header,RNX_data]=load_rinex_o(rinex_o);
if RNX_header.VERSION==0
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

% add SA debug
%~ [RNX_header,RNX_data] = tool_add_SA(RNX_header,RNX_data);

%%%%% Set epoch min and epoch max
% Epoch_min
if(Epoch_min < 1 || Epoch_min > size(RNX_data.G,3))
	Epoch_min = 1;
end

% Number of epochs
if(N_epoch < 1 || N_epoch > size(RNX_data.G,3)-Epoch_min+1)
	N_epoch = size(RNX_data.G,3)-Epoch_min+1;
end


%%%%% Initial values
% X0 if not specified in options
if(sum(X0)==0)
	X0 = [RNX_header.X;RNX_header.Y;RNX_header.Z]; % initial position
end  

% station approx position
sta_pos = [X0(1:3)];

% phase center and antenna height correction -> approx pos == phase center
% calibrations
antex_file = 'igs08.atx';
[ATX_header,ATX_data] = load_antex(antex_file);
[ATX1] = get_antex(ATX_header, ATX_data, RNX_header.ANT_TYPE, 'G01');
[ATX2] = get_antex(ATX_header, ATX_data, RNX_header.ANT_TYPE, 'G02');

if strcmp(freq, 'F1')            % First frequency
	dE = 0;
	dN = 0;
	dU = ATX1.UP;
elseif strcmp(freq, 'F2')        % Second frequency
	dE = 0;
	dN = 0;
	dU = ATX2.UP;
else                             % Iono-free
	eph_atx.const = 'G';
	dE = 0;
	dN = 0;
	dU = 2.5*ATX1.UP-1.5*ATX2.UP;
end

% antenna height
dE = dE + RNX_header.dE;
dN = dN + RNX_header.dN;
dU = dU + RNX_header.dH;

% station position to phase center (approximate coordinates)
[X0(1),X0(2),X0(3)] = corr_pos_atx(X0(1),X0(2),X0(3), dE, dN, dU, 1);

% initial cdtr
X0 = [X0(1:3);0];                                       

% initial cGGTO and cGPGL
if (strcmp(nav,'brdc'))                            

	if isfield(NAV_header,'GPGA')
		cGGTO = c * NAV_header.GPGA(1);  % Approx GGTO, not the real formula which need mjd (cf official doc Galileo)
		X0 = [X0;cGGTO];  
	else
		X0 = [X0;0];  
	end

	if isfield(NAV_header,'GLGP')
		cGPGL = c * NAV_header.GLGP(1); % in rinex, GPGL is not defined, but  - GLGP is present
		X0 = [X0;cGPGL];  
	else
		X0 = [X0;0];  
	end
	
else % sp3
	X0 = [X0;0;0]; 
end

%%%%% Set options
options = cell(0);

comp_options.X0 = X0;
comp_options.freq = freq;
comp_options.iono = iono;
comp_options.nav = nav;
comp_options.const = const;
comp_options.sat_num = sat_num;
comp_options.verbose = verbose;
comp_options.cut_off = cut_off;
comp_options.global_pos = global_pos;
comp_options.Epoch_min = Epoch_min;
comp_options.N_epoch = N_epoch;


%%%%% Displaying initial values
tool_print_info('',1);
tool_print_info('INITIAL VALUES : ',1);
tool_print_info(sprintf('\tINITIAL COORDINATES (M): %0.3f %0.3f %0.3f',X0(1:3)),1);
tool_print_info(sprintf('\tINITIAL cDTR (M) : %0.3f', X0(4)),1);
tool_print_info(sprintf('\tINITIAL cGGTO (M) : %0.3f', X0(5)),1);
tool_print_info(sprintf('\tINITIAL cGPGL (M) : %0.3f', X0(6)),1);
tool_print_info('',1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute one pos, dtr, GGTO and GPGL per epoch
if global_pos==0

	global_process_time = now;
	mean_t = 0;
	
	for i = 1:N_epoch
			
		tool_print_info('',0)
		tool_print_info(sprintf('EPOCH %d/%d : ',i,N_epoch),1)
		tool_print_info(sprintf('\tEST REM TIME : %2d MIN %0d S   ',floor((N_epoch-i)*mean_t/60),fix(((N_epoch-i)*mean_t/60-floor((N_epoch-i)*mean_t/60))*60)),0);
	
		process_time = now;
						
		% Preprocessing
		[G,calc]=calc_preprocessing_spp(RNX_header,RNX_data,NAV_header,NAV_data,Epoch_min+i-1,comp_options); % begin at Epoch_min
		t = calc.t.mjd*ones(calc.nb_sat,1);
		PosSat = calc.PosSat;
		Dobs = calc.Dobs;
		ElevSat = calc.ElevSat;
		sat_index = calc.sat_index;
		nb_loc_val = calc.nb_sat;
						
		tool_print_info(sprintf('\tSAT : %2d GPS, %2d GLO, %2d GAL',calc.nb_GPS,calc.nb_GLO,calc.nb_GAL),0);
		
		% Least Square computation	
		[result_i] = calc_LS_code(t,PosSat,Dobs,ElevSat,sat_index,X0);
		
		% Statistic indicators computation
		[stat]=calc_stat_indic(result_i.Qxx,result_i.X,result_i.Y,result_i.Z);
		
		process_time = (now - process_time)*86400;
		
		mean_t = mean([process_time mean_t]);
		
		% Update initial coordinates + cdtr
		X0(1:4) = [result_i.X;result_i.Y;result_i.Z;result_i.cdtr];
		options.X0 = X0;
		
		% phase center and antenna height correction
		
		[result_i.X, result_i.Y, result_i.Z] = corr_pos_atx(result_i.X, result_i.Y, result_i.Z, dE, dN, dU,2);

	
		% Save results and intermediates values -> for matlab, need to preallocate result 
		% otherwise error 'the input character is not valid in MATLAB statements or expressions'
		if i == 1 
			interm_results.G = G;
			interm_results.calc = calc;
			result = tool_save_res(sta_pos,result_i,stat,process_time);	

		else
			interm_results(i).G = G;
			interm_results(i).calc = calc;
			result(i) = tool_save_res(sta_pos,result_i,stat,process_time);	

		end
		
				
	end
	
	tool_print_info('',1);
	
	%%%%% plots
		
	% plot results
	plot_results_code(result,rinex_o(end-11:end-8),dir_out);

	% skyplot
	if ~strcmp(dir_out,'')
		save_skyplot = strcat(dir_out,filesep(),'skyplot');
	else
		save_skyplot = '';
	end	
	[mat_skyplot,sat_index] = tool_prep_skyplot(interm_results);	
	plot_skyplot(mat_skyplot,sat_index,result(1).sta_pos(1),result(1).sta_pos(2),result(1).sta_pos(3),rinex_o(end-11:end-8),save_skyplot);
	
	global_process_time = (now - global_process_time)*86400;
	
	tool_print_info(' ',1);

	% report
	if ~strcmp(dir_out,'')
		save_report = strcat(dir_out,filesep(),'results_',datestr(now,'yy_mm_dd_HHMMss'),'.txt');
	else
		save_report = strcat('results_',datestr(now,'yy_mm_dd_HHMMss'),'.txt');
	end
	tool_save_report(result,save_report);
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute one pos, GGTO and GPGL for all rinex file, and one dtr per epoch
elseif global_pos==1

	global_process_time = now;
	
	% Preprocessing : satellite position and corrected PR computation
		
	tool_print_info(' ',1);
	tool_print_info('PREPROCESSING',1);
	tool_print_info('',1);
		
	
	% matrix initialisation
	PosSat = zeros(sat_num*N_epoch,3);
	Dobs = zeros(sat_num*N_epoch,1);
	AzSat = zeros(sat_num*N_epoch,1);
	ElevSat = zeros(sat_num*N_epoch,1);
	sat_index = cell(sat_num*N_epoch,2);
	t = zeros(sat_num*N_epoch,1);
	N_val = 0;
	
	
	for i = 1:N_epoch
	
		tool_print_info('',0)
		tool_print_info(sprintf('EPOCH %d/%d : ',i,N_epoch),1)

		% Preprocessing
		[G,calc]=calc_preprocessing_spp(RNX_header,RNX_data,NAV_header,NAV_data,Epoch_min+i-1,comp_options); % begin at epoch_min
		
		tool_print_info(sprintf('\tSAT : %2d GPS, %2d GLO, %2d GAL',calc.nb_GPS,calc.nb_GLO,calc.nb_GAL),0);
		
		nb_loc_val = size(calc.PosSat,1);

		PosSat(N_val+1:N_val + nb_loc_val,:) = calc.PosSat;
		Dobs(N_val+1:N_val + nb_loc_val,:) = calc.Dobs;
		AzSat(N_val+1:N_val + nb_loc_val,:) = calc.AzSat;
		ElevSat(N_val+1:N_val + nb_loc_val,:) = calc.ElevSat;
		
		sat_index(N_val+1:N_val + nb_loc_val,1) = calc.sat_index(:,1);
		sat_index(N_val+1:N_val + nb_loc_val,2) = calc.sat_index(:,2);
			
		t(N_val+1:N_val + nb_loc_val,:) = calc.t.mjd*ones(nb_loc_val,1);

		N_val = N_val + nb_loc_val;
					
		interm_results(i).G = G;
		interm_results(i).calc = calc;
			
	end

	% redim matrix
	PosSat = PosSat(1:N_val,:);
	Dobs = Dobs(1:N_val,:);
	AzSat = AzSat(1:N_val,:);
	ElevSat = ElevSat(1:N_val,:);
	sat_index = sat_index(1:N_val,:);
	t = t(1:N_val,:);
				
	tool_print_info('',1);
	tool_print_info('LEAST SQUARE COMPUTATION',1);
	tool_print_info('',1);
				
	% Least Square computation	
	[result_i] = calc_LS_code(t,PosSat,Dobs,ElevSat,sat_index,X0);
		
	% Statistic indicators computation
	[stat]=calc_stat_indic(result_i.Qxx,result_i.X,result_i.Y,result_i.Z);
	
	% phase center and antenna height correction
	[result_i.X, result_i.Y, result_i.Z] = corr_pos_atx(result_i.X, result_i.Y, result_i.Z, dE, dN, dU,2);

	% save results
	[result] = tool_save_res(sta_pos,result_i,stat,global_process_time);
	
	tool_print_info('',1);

	%%%%% plot
	if ~strcmp(dir_out,'')
		save_hist_res = strcat(dir_out,filesep(),'hist_res');
	else
		save_hist_res = '';
	end
	figure()
	plot_hist(result.Vnorm,40,'Histogram : normalized residuals','Value',save_hist_res);

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
		save_skyplot = strcat(dir_out,filesep(),'skyplot');
	else
		save_skyplot = '';
	end	
	mat_skyplot = [t,AzSat,ElevSat];
	plot_skyplot(mat_skyplot,sat_index,result.X,result.Y,result.Z,rinex_o(end-11:end-8),save_skyplot);
	
	global_process_time = (now - global_process_time)*86400;
	
	% report
	if ~strcmp(dir_out,'')
		save_report = strcat(dir_out,filesep(),'results_',datestr(now,'yy_mm_dd_HHMMss'),'.txt');
	else
		save_report = strcat('results_',datestr(now,'yy_mm_dd_HHMMss'),'.txt');
	end
	tool_save_report(result,save_report);
	
end

tool_print_info('',1);
tool_print_info(sprintf('TOTAL PROCESSING TIME : %0.3f S',global_process_time),1);	
tool_print_info('',1);


tool_print_info('----------------------------------------',1);
tool_print_info('         END FUNCTION RUN_SPP',1);
tool_print_info('----------------------------------------',1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [options,X0,freq,cut_off,iono,verbose,nav,global_pos,Epoch_min,N_epoch,const,sat_num,dir_out] = set_options(options,X0,freq,cut_off,iono,verbose,nav,global_pos,Epoch_min,N_epoch,const,sat_num,dir_out)
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
    
    if isfield(options,'global_pos') 
        if (options.global_pos==1)
            global_pos = 1;
        else 
			options.global_pos = global_pos;
        end  
    else 
		options.global_pos = global_pos;
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
