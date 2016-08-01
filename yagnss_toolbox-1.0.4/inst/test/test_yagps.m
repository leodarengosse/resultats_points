function test_yagps(test_num)
%% function yagps(num_test)
%% Run tests
%%
%% Clement Fontaine 2013-12-02
%% 
%% Input : 
%% - test_num : test number : 
%%
%% Test 1 : spp smne, G, brdc, iono_free
%% Test 2 : spp mlvl, G, brdc, iono_free, global   
%% Test 3 : spp smne, GE, brdc, iono_free
%% Test 4 : spp smne, GRE, brdc, iono_free
%% Test 5 : spp smne, G, sp3, iono_free
%% Test 6 : spp smne, G, brdc, F1
%% Test 7 : spp smne, G, brdc, klobuchar
%% Test 8 : spp smne, E, brdc, iono_free
%% Test 9 : spp tlse, GE, brdc, iono_free
%% Test 10 : DGNSS smne->mlvl, GE, brdc, F1       
%% Test 11 : DGNSS smne->mlvl, GE, brdc, iono_free       
%% Test 12 : DGNSS smne->mlvl, GE, sp3, F1       
%% Test 13 : DGNSS tlmf->tlse, GE, brdc, F1       
%% Test 14 : DGNSS smne->mlvl, E, brdc, F1       
%% Test 15 : phase smne->mlvl G
%% Test 16 : phase smne->mlvl GE, 100
%% Test 17 : phase smne->mlvl GRE, 100
%% Test 18 : phase smne->mlvl E, 100
%% Test 19 : phase smne->mlvl GRE, 100
%% Test 20 : phase smne->mlvl G sp3
%% Test 21 : phase smne->mlvl GRE
%% Test 22 : phase tlmf->tlse GE, brdc
%% Test 23 : plot_orbit
%% Test 24 : comp sp3-brdc for sat G22
%% Test 25 : comp igr-igs for sat G22
%% Test 26 : Comp sp3
%% Test 27 : ANTEX   
%% Test 28 : run_plot_obs  
%%
%% 'test_yagps' without args runs all the tests  
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
more off
%% format long

global VERBOSE;
VERBOSE = 1;

brdm = 'brdm1500.13p';
%~ brdm = 'brdm1300.13p';
sp3 = {'igs17424.sp3','igl17424.sp3','grm17424.sp3'};
igr = {'igr17424.sp3'};
%~ sp3 = {'igs17395.sp3','igl17395.sp3','grm17395.sp3'};
atx = 'igs08.atx';

mjd_min = 56442.0;
mjd_max = 56443.0;
%~ mjd_min = 56422.0;
%~ mjd_max = 56423.0;

rinex_mlvl = 'mlvl1500.13o';
rinex_smne = 'smne1500.13o';
rinex_tlse = 'tlse1500.13o';
rinex_tlmf = 'tlmf1500.13o';
%~ rinex_mlvl = 'mlvl1300.13o';
%~ rinex_smne = 'smne1300.13o';
%~ rinex_tlse = 'tlse1300.13o';
%~ rinex_tlmf = 'tlmf1300.13o';
X0_smne = [4201792.2950;177945.2380;4779286.6850]; 
X0_mlvl = [4201577.2090;189859.8560;4779064.5670];
X0_tlse = [4627852.0690;119639.7490;4372993.3260];
X0_tlmf = [4627010.0560;111069.2130;4374122.8590];

if nargin == 1
	switch test_num 
	
		case 1
			tool_print_info('Test 1 : spp smne, G, brdc, iono_free',1);
				
			options.verbose = 1;
			options.const = 'G';
			options.nav='brdm';
			options.freq = 'iono_free';
			options.cut_off=7;
			options.global_pos=0;
			options.Epoch_min = 1;
			options.N_epoch = 800;
			options.dir_out = 'Test_1_spp_smne_G_brdc_ionofree';
			
			run_spp(rinex_smne,brdm,options);
			
		case 2
	        tool_print_info('Test 2 : spp mlvl, G, brdc, iono_free, global',1);       
	        
	   		options.verbose = 1;
			options.const = 'G';
			options.nav='brdm';
			options.freq = 'iono_free';
			options.cut_off=7;
			options.global_pos=1;
			options.Epoch_min = 1;
			options.N_epoch = 50;
			options.dir_out = 'Test_2_spp_smne_G_brdc_ionofree_global';
			
			run_spp(rinex_smne,brdm,options);
			
		case 3
	        tool_print_info('Test 3 : spp smne, GE, brdc, iono_free',1);
	        
	   		options.verbose = 1;
			options.const = 'GE';
			options.nav='brdm';
			options.freq = 'iono_free';
			options.cut_off=7;
			options.global_pos=0;
			options.Epoch_min = 1;
			options.N_epoch = 800;
			options.dir_out = 'Test_3_spp_smne_GE_brdc_ionofree';
			
			run_spp(rinex_smne,brdm,options);
			
		case 4
	        tool_print_info('Test 4 : spp smne, GRE, brdc, iono_free',1);
	        
	   		options.verbose = 1;
			options.const = 'GRE';
			options.nav='brdm';
			options.freq = 'iono_free';
			options.cut_off=7;
			options.global_pos=0;
			options.Epoch_min = 1;
			options.N_epoch = 800;
			options.dir_out = 'Test_4_spp_smne_GRE_brdc_ionofree';
			
			run_spp(rinex_smne,brdm,options);
			
		case 5
	        tool_print_info('Test 5 : spp smne, G, sp3, iono_free',1);
	        
	   		options.verbose = 1;
			options.const = 'G';
			options.nav='sp3';
			options.freq = 'iono_free';
			options.cut_off=7;
			options.global_pos=0;
			options.Epoch_min = 1;
			options.N_epoch = 800;
			options.dir_out = 'Test_5_spp_smne_G_sp3_ionofree';
			
			run_spp(rinex_smne,sp3,options);
			
			
	    case 6
	        tool_print_info('Test 6 : spp smne, G, brdc, F1',1);
	        
	   		options.verbose = 1;
			options.const = 'G';
			options.nav='brdm';
			options.freq = 'F1';
			options.cut_off=7;
			options.global_pos=0;
			options.Epoch_min = 1;
			options.N_epoch = 800;
			options.dir_out = 'Test_6_spp_smne_G_brdc_F1';
			
			run_spp(rinex_smne,brdm,options);
								
			
	    case 7
	        tool_print_info('Test 7 : spp smne, G, brdc, klobuchar',1);
	        
	   		options.verbose = 1;
			options.const = 'G';
			options.nav='brdm';
			options.freq = 'F1';
			options.iono = 'klobuchar';
			options.cut_off=7;
			options.global_pos=0;
			options.Epoch_min = 1;
			options.N_epoch = 800;
			options.dir_out = 'Test_7_spp_smne_G_brdc_klobuchar';
			
			run_spp(rinex_smne,brdm,options);
			
	    case 8
	        tool_print_info('Test 8 : spp smne, E, brdc, iono_free',1);
	        
	   		options.verbose = 1;
			options.const = 'E';
			options.nav='brdm';
			options.freq = 'iono_free';
			options.cut_off=7;
			options.global_pos=0;
			options.Epoch_min = 1;
			options.N_epoch = 200;
			options.dir_out = 'Test_8_spp_smne_E_brdc_iono_free';
			
			run_spp(rinex_smne,brdm,options);
			
	    case 9
	    
	        tool_print_info('Test 9 : spp tlse, GE, brdc, iono_free',1);
	        
	   		options.verbose = 1;
			options.const = 'GE';
			options.nav='brdm';
			options.freq = 'iono_free';
			options.cut_off=7;
			options.global_pos=0;
			options.Epoch_min = 1;
			options.N_epoch = 200;
			options.dir_out = 'Test_9_spp_tlse_GE_brdc_iono_free';
			
			run_spp(rinex_tlse,brdm,options);
	

		case 10
	        tool_print_info('Test 10 : DGNSS smne->mlvl, GE, brdc, F1',1);       
	        
			options.verbose = 1;
			options.const = 'GE';
			options.nav='brdc';
			options.freq='F1';
			options.cut_off=7;
			options.Epoch_min = 1;
			options.N_epoch = 800;
			options.dir_out = 'Test_10_DGNSS_smne_mlvl_GE_brdc_F1';
				
			run_DGNSS(rinex_smne,rinex_mlvl,brdm,X0_smne,options);
			
	    case 11
	        tool_print_info('Test 11 : DGNSS smne->mlvl, GE, brdc, iono_free',1);       
	        
			options.verbose = 1;
			options.const = 'GE';
			options.nav='brdc';
			options.freq='iono_free';
			options.cut_off=7;
			options.Epoch_min = 1;
			options.N_epoch = 800;
			options.dir_out = 'Test_11_DGNSS_smne_mlvl_GE_brdc_ionofree';
				
			run_DGNSS(rinex_smne,rinex_mlvl,brdm,X0_smne,options);
			
	    case 12
	        tool_print_info('Test 12 : DGNSS smne->mlvl, GE, sp3, F1',1);       
	        
			options.verbose = 1;
			options.const = 'GE';
			options.nav='sp3';
			options.freq='F1';
			options.cut_off=7;
			options.Epoch_min = 1;
			options.N_epoch = 800;
			options.dir_out = 'Test_12_DGNSS_smne_mlvl_GE_sp3_ionofree';
				
			run_DGNSS(rinex_smne,rinex_mlvl,sp3,X0_smne,options);
			
			
		case 13
	        tool_print_info('Test 13 : DGNSS tlmf->tlse, GE, brdc, F1',1);       
	        
			options.verbose = 1;
			options.const = 'GE';
			options.nav='brdc';
			options.freq='F1';
			options.cut_off=7;
			options.Epoch_min = 1;
			options.N_epoch = 800;
			options.dir_out = 'Test_13_DGNSS_tlmf_tlse_GE_brdc_F1';
				
			run_DGNSS(rinex_tlmf,rinex_tlse,brdm,X0_tlmf,options);
			
		case 14
	        tool_print_info('Test 14 : DGNSS smne->mlvl, E, brdc, F1',1);       
	        
			options.verbose = 1;
			options.const = 'E';
			options.nav='brdc';
			options.freq='F1';
			options.cut_off=7;
			options.Epoch_min = 1;
			options.N_epoch = 200;
			options.dir_out = 'Test_14_DGNSS_smne_mlvl_E_brdc_F1';
				
			run_DGNSS(rinex_smne,rinex_mlvl,brdm,X0_smne,options);
			
			
		case 15
			tool_print_info('Test 15 : phase smne->mlvl G',1);
			
			options.const = 'G';
			options.cut_off = 7;
			options.dir_out = 'Test_15_phase_smne_mlvl_brdc_G';
			options.Epoch_min = 1;
			options.N_epoch = 800;
			
			run_phase(rinex_smne,rinex_mlvl,brdm,X0_smne,options);
			
		case 15
			tool_print_info('Test 15 : phase smne->mlvl G, 100',1);
			
			options.const = 'G';
			options.cut_off = 7;
			options.dir_out = 'Test_15_phase_smne_mlvl_brdc_G_100';
			options.Epoch_min = 1;
			options.N_epoch = 100;
			
			run_phase(rinex_smne,rinex_mlvl,brdm,X0_smne,options);			
			
	    case 16
			tool_print_info('Test 16 : phase smne->mlvl GE, 100',1);
			
			options.const = 'GE';
			options.cut_off = 7;
			options.dir_out = 'Test_16_phase_smne_mlvl_GE_brdc_100';
			options.Epoch_min = 1;
			options.N_epoch = 100;
			
			run_phase(rinex_smne,rinex_mlvl,brdm,X0_smne,options);		
				
	    case 17
			tool_print_info('Test 17 : phase smne->mlvl GRE, 100',1);
			
			options.const = 'GRE';
			options.cut_off = 7;
			options.dir_out = 'Test_17_phase_smne_mlvl_GRE_brdc_100';
			options.Epoch_min = 1;
			options.N_epoch = 100;
			
			run_phase(rinex_smne,rinex_mlvl,brdm,X0_smne,options);			
			
	    case 18
			tool_print_info('Test 18 : phase smne->mlvl E, 100',1);
			
			options.const = 'E';
			options.cut_off = 7;
			options.dir_out = 'Test_18_phase_smne_mlvl_E_brdc_100';
			options.Epoch_min = 1;
			options.N_epoch = 100;
			
			run_phase(rinex_smne,rinex_mlvl,brdm,X0_smne,options);	
		
		case 19
			tool_print_info('Test 19 : phase smne->mlvl GRE, 100',1);
			
			options.const = 'GRE';
			options.cut_off = 7;
			options.dir_out = 'Test_19_phase_smne_mlvl_GRE_brdc_100';
			options.Epoch_min = 1;
			options.N_epoch = 100;
			
			run_phase(rinex_smne,rinex_mlvl,brdm,X0_smne,options);	
			
	
		case 20
			tool_print_info('Test 20 : phase smne->mlvl G sp3',1);
			
			options.nav = 'sp3';
			options.const = 'G';
			options.cut_off = 7;
			options.dir_out = 'Test_20_phase_smne_mlvl_G_sp3';
			options.Epoch_min = 1;
			options.N_epoch = 250;
			
			run_phase(rinex_smne,rinex_mlvl,sp3,X0_smne,options);
			
		case 21
			tool_print_info('Test 21 : phase smne->mlvl GRE',1);
			
			options.const = 'GRE';
			options.cut_off = 7;
			options.dir_out = 'Test_21_phase_smne_mlvl_GRE_brdc';
			options.Epoch_min = 1;
			options.N_epoch = 800;
			
			run_phase(rinex_smne,rinex_mlvl,brdm,X0_smne,options);	

	    case 22
			tool_print_info('Test 22 : phase tlmf->tlse GE, brdc',1);
			
			options.const = 'GE';
			options.cut_off = 7;
			options.dir_out = 'Test_22_phase_tlmf_tlse_GE_brdc';
			options.Epoch_min = 1;
			options.N_epoch = 800;
			
			run_phase(rinex_tlmf,rinex_tlse,brdm,X0_tlmf,options);		
			
		case 23
			tool_print_info('Test 23 :  plot_orbit',1)
		
			savef = 'Test_23_plot_orbit';
			sat = 'G19G21G22G23G24G25';
			[NAV_header,NAV_data] = load_rinex_n(brdm);
			
			plot_orbit(NAV_header,NAV_data,sat,mjd_min,mjd_max,savef);			
	         
	    case 24
	        tool_print_info('Test 24 : comp sp3-brdc for sat G22',1);
	        
			step = 30;
			order = 9;
			savef = 'Test_24_comp_sp3_brdc_G22';
			sat = 'G22';
	
			run_comp_nav_sp3(brdm,sp3,sat,mjd_min,mjd_max,step,order,savef)
			
	    case 25
	        tool_print_info('Test 25 : comp igr igs sat G22',1);
	        
			step = 30;
			order = 9;
			savef = 'Test_25_comp_igr_igs_G22';
			sat = 'G22';
	
			run_comp_2_sp3(igr,sp3,sat,mjd_min,mjd_max,step,order,savef);
			
			
		case 26
	        tool_print_info('Test 26 : Comp sp3',1);
	        
	        run_comp_sp3(sp3, 'G14',mjd_min+0.1, mjd_max-0.1,120,'Test_25_comp_sp3_G14')
	
			
		case 27
	        tool_print_info('Test 27 : ANTEX',1);       
	         
	        savef = 'Test_26_plot_antex';
	        [ATX_header,ATX_data]=load_antex(atx);
	        plot_antex(ATX_header,ATX_data,'TIAPENG3100R1   NONE','G01',savef);
	    
	    case 28
			tool_print_info('Test 28 : run_plot_obs',1);
			
			 run_plot_obs(rinex_smne, 'E', 12, 'Test 28 : obs_sat');
	        		
	    otherwise
	        tool_print_info('Run all',1);
		    for i = 1:28
				test_yagps(i);
		    end
	        toc
	end

	else
			tic
		    tool_print_info('Run all',1);
		    for i = 1:28
				test_yagps(i);
		    end
	        toc
	
	end



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
