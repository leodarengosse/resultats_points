function [Eph]=get_ephemeris(NAV_header,NAV_data,constellation,PRN,mjd)
%% function [Eph]=get_ephemeris(NAV_header,NAV_data,constellation,PRN,mjd)
%% Get navigation message for one satellite at mjd (GPS, Glonass and Galileo supported)
%%
%% Jacques Beilin - ENSG/DPTS - 2010-02-25
%% Clement Fontaine 2013-10-16
%%
%% Input :
%% - NAV_header : structure containing navigation Rinex header
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
%% - NAV_data : matrix containing data
%%   NAV_header and NAV_data are set up with function load_rinex_n
%% - constellation : 'G' = GPS, 'R' = Glonass, 'E' = Galileo
%% - PRN : satellite id
%% - mjd : date (modified julian date format)
%%
%% Output : 
%% 	- Eph structure containing informations
%%  
%% If no informations are found, no fields are defined in Eph. 
%% Ex : isfield(Eph,'mjd') returns 0 if Eph is empty.
%%
%% Eph content depends on constellation (GPS, Glonass and Galileo supported)
%%
%% Note : get_ephemeris returns 0 if health of satellite is not OK (sv_health ~= 0)
%% 
%%
%% Examples :
%% 
%% GPS
%%
%% [Eph]=get_ephemeris(NAV_header,NAV_data,'G',1,56442.0833333)
%% Eph =
%% {
%%   PRN =  1
%%   mjd =  56442.0833333335
%%   TOC =  56442.0833333335
%%   alpha0 =  3.15997749567000e-05
%%   alpha1 =  4.32009983342100e-12
%%   alpha2 = 0
%%   const = G
%%   IODE =  94
%%   crs =  2.25000000000000
%%   delta_n =  4.63876465173900e-09
%%   M0 =  1.81168558621700
%%   cuc =  2.75671482086200e-07
%%   e =  0.00189499533735200
%%   cus =  1.05910003185300e-05
%%   sqrt_a =  5153.70066070600
%%   TOE =  352800
%%   cic =  2.60770320892300e-08
%%   OMEGA = -1.46537454108600
%%   cis =  1.37835741043100e-07
%%   i0 =  0.959848465747100
%%   crc =  173.687500000000
%%   omega =  0.255929757187400
%%   OMEGA_DOT = -8.04104922769100e-09
%%   IDOT =  3.31799535068200e-10
%%   code_L2 = 0
%%   gps_wk =  1742
%%   L2_P = 0
%%   sv_acc =  2
%%   sv_health = 0
%%   TGD =  8.38190317153900e-09
%%   IODC =  94
%%   trans_time =  345600
%% }
%%
%%  GLONASS
%%
%% [Eph]=get_ephemeris(NAV_header,NAV_data,'R',1,56442.0833333)
%% Eph =
%% {
%%   PRN =  1
%%   mjd =  56442.0104166665
%%   TOC =  56442.0104166665
%%   SV_clock_offset = -1.72349624335800e-04
%%   SV_relat_freq_offset = 0
%%   const = R
%%   Message_frame_time =  345600
%%   X =  18585.5019531200
%%   X_dot =  1.61786460876500
%%   X_acc = 0
%%   sv_health = 0
%%   Y = -12058.0571289100
%%   Y_dot = -0.623869895935100
%%   Y_acc =  9.31322574615500e-10
%%   freq_num =  1
%%   Z = -12625.5478515600
%%   Z_dot =  2.97767543792700
%%   Z_acc = 0
%%   age_op_inf = 0
%% }
%% 
%%  Galileo
%% 
%% [Eph]=get_ephemeris(NAV_header,NAV_data,'E',11,56442.0833333)
%% Eph =
%% {
%%   PRN =  11
%%   mjd =  56442
%%   TOC =  56442
%%   alpha0 =  9.06531931832400e-04
%%   alpha1 =  8.25934876047500e-11
%%   alpha2 = 0
%%   const = E
%%   IODnav =  64
%%   crs = -127
%%   delta_n =  3.22049128924600e-09
%%   M0 =  2.45252068075900
%%   cuc = -5.72949647903400e-06
%%   e =  3.02415806800100e-04
%%   cus =  1.07511878013600e-05
%%   sqrt_a =  5440.61737251300
%%   TOE =  345600
%%   cic =  2.42143869400000e-08
%%   OMEGA = -2.27501476678800
%%   cis = 0
%%   i0 =  0.957955133762200
%%   crc =  109.468750000000
%%   omega = -0.819738964252500
%%   OMEGA_DOT = -5.59094717110500e-09
%%   IDOT = -1.00004165574900e-11
%%   data_src =  513
%%   gal_wk =  1742
%%   SISA = -1
%%   sv_health =  452
%%   BGDE5a = -6.51925802230800e-09
%%   BGDE5b = 0
%%   trans_time =  346255
%% }
%% 
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Eph=[];
%Eph.PRN = 0;
delta_t=2/24; % max delta of 2 hours
min_val = 10;


if strcmp(constellation,'G')
	NAV_data = NAV_data.G;
elseif strcmp(constellation,'R')
	NAV_data = NAV_data.R;
elseif strcmp(constellation,'E')
	NAV_data = NAV_data.E;
else
	tool_print_info('Constellation not implemented : Eph = []',2);
	return
end

if (PRN>32 || PRN <1)
	tool_print_info('PRN < 1 or PRN > 32 : Eph = []',2);
	return;
end

[Nsat,Nfield,Nepoch]=size(NAV_data);

% 1 - Select Eph with Eph.mjd close to mjd

% min_val -> index of closest ephemeris in NAV_data
min_val = find(squeeze(abs(mjd - NAV_data(PRN,1,:)))==min(squeeze(abs(mjd - NAV_data(PRN,1,:)))));
data = squeeze(NAV_data(PRN,:,min_val));

if(min(size(data))~=1)
	data = data';
end

% no eph found
if length(min_val)==0
	return;
end

% test : mjd close to TOE
if(abs(mjd - NAV_data(PRN,1,min_val(1)))>=delta_t)
	tool_print_info('TOE is too far from mjd : Eph = []',4);
	return;
end

% if several Eph, select Eph for which sv_health is Ok (=0)
if strcmp(constellation,'G')

	i = find(data(:,26) == min(data(:,26)));
	i = i(1);
	
	if ~data(i,26)== 0 % health not ok
		return;
	else
		data = data(i,:);
	end
	
elseif strcmp(constellation,'E')


	% 2 - If 3 mjd, select latest date
	
	mjd_max = max(data(:,1));
	ind_data = find(data(:,1) == mjd_max);
	if ind_data~=0
		data = data(ind_data,:);
	end
		
	% select data with health code ok	
	
	health = data(:,25);
	src = data(:,22);
	[num] = test_health_GAL(health,src);
	data = data(num,:);
	

	% median of remaining ephemeris
	
	if size(data,1)==0
		return;
	end
	
	data = median(data,1);

	
	
	
elseif strcmp(constellation,'R')

	i = find(data(:,8) == min(data(:,8)));
	i = i(1);
	
	if ~data(i,8)== 0 % health not ok
		return;
	else
		data = data(i,:);
	end
	

	
end

		
Eph.PRN 	  = PRN;
Eph.mjd       = data(1);
Eph.TOC       = Eph.mjd;

Eph.alpha0    = data(2);
Eph.alpha1    = data(3);
Eph.alpha2    = data(4);

if Nfield==29 % GPS navigation message

	Eph.const = 'G';
	Eph.IODE      = data(5);
	Eph.crs       = data(6);
	Eph.delta_n   = data(7);
	Eph.M0        = data(8);

	Eph.cuc       = data(9);
	Eph.e         = data(10);
	Eph.cus       = data(11);
	Eph.sqrt_a    = data(12);

	Eph.TOE       = data(13);
	Eph.cic       = data(14);
	Eph.OMEGA     = data(15);
	Eph.cis       = data(16);

	Eph.i0        = data(17);
	Eph.crc       = data(18);
	Eph.omega     = data(19);
	Eph.OMEGA_DOT = data(20);

	Eph.IDOT      = data(21);
	Eph.code_L2   = data(22);
	Eph.gps_wk    = data(23);
	Eph.L2_P      = data(24);

	Eph.sv_acc    = data(25);
	Eph.sv_health = data(26);
	Eph.TGD       = data(27);
	Eph.IODC      = data(28);

	Eph.trans_time = data(29);
	

elseif Nfield == 28 % Galileo navigation message

	Eph.const = 'E';
	Eph.IODnav    = data(5);
	Eph.crs       = data(6);
	Eph.delta_n   = data(7);
	Eph.M0        = data(8);

	Eph.cuc       = data(9);
	Eph.e         = data(10);
	Eph.cus       = data(11);
	Eph.sqrt_a    = data(12);
                    
	Eph.TOE       = data(13);
	Eph.cic       = data(14);
	Eph.OMEGA     = data(15);
	Eph.cis       = data(16);
                    
	Eph.i0        = data(17);
	Eph.crc       = data(18);
	Eph.omega     = data(19);
	Eph.OMEGA_DOT = data(20);

	Eph.IDOT      = data(21);
	Eph.data_src  = data(22);
	Eph.gal_wk    = data(23);

	Eph.SISA      = data(24);
	Eph.sv_health = data(25);
	Eph.BGDE5a    = data(26);
	Eph.BGDE5b    = data(27);

	Eph.trans_time = data(28);	

elseif Nfield == 16 % GLONASS navigation message
	
	Eph.const='R';
	Eph = rmfield(Eph,'alpha2');
	Eph.SV_clock_offset = Eph.alpha0;
	Eph = rmfield(Eph,'alpha0');
	Eph.SV_relat_freq_offset = Eph.alpha1;
	Eph = rmfield(Eph,'alpha1');


	Eph.Message_frame_time   = data(4);;
	
	Eph.X          = data(5)*10^3;
	Eph.X_dot      = data(6)*10^3;
	Eph.MS_X_acc   = data(7)*10^3;
	Eph.sv_health     = data(8);
	               
	Eph.Y          = data(9)*10^3;
	Eph.Y_dot      = data(10)*10^3;
	Eph.MS_Y_acc   = data(11)*10^3;
	Eph.freq_num   = data(12);
	                 
	Eph.Z          = data(13)*10^3;
	Eph.Z_dot      = data(14)*10^3;
	Eph.MS_Z_acc   = data(15)*10^3;
	Eph.age_op_inf = data(16);
   
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [num] = test_health_GAL(health,src)
%% function [num] = test_health_GAL(health)
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%~ exp = '^0(00|11)0(00|11)...$';
%~ exp = '^0(00|11).(00|11)...$';
%~ exp = '^000000...$';
%~ exp = '(000)(000|111)...';
%~ exp = '000000000';
%~ exp = '000000...';
%~ exp = '000(000|111)...';
%~ exp = '(000|111)000...';
%~ exp = '(000|111)111...';


exp = '^(0|1)(11|00)(0|1)(11|00)...$';

% dec 2 bin

num = [];

for i = 1:length(health)




	%~ num = i;
	%~ return;
	
	if health(i) == 0
		num = i;
		return;
	end
	
	bin_h = dec2bin(health(i),9);
	

	match = regexp(bin_h,exp);
	if (match==1)
		num = [num;i];
	end

	

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

