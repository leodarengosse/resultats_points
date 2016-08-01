function [X, Y, Z, DDN, V, Vnorm, sigma02, Qxx] = interf_calc_LS_phase(t,Pobs,PosSat,ElevSat,Dtropo,sat_const, sat_PRN, amb_index,amb0,X0_base,X0_rover,pivots,interv_pivot)
%% function [X, Y, Z, DDN, V, Vnorm, sigma02, Qxx] = interf_calc_LS_phase(t,Pobs,PosSat,ElevSat,Dtropo,sat_const, sat_PRN, amb_index,amb0,X0_base,X0_rover,pivots,interv_pivot)
%%
%% Rover position and ambiguity estimation for phase processing -> interface function
%% P depends on elevation of satellites
%%
%% Clement Fontaine 2014-01-29
%%
%% Input : 
%% - t : vector containing time of obs (one time per obs) [t] (mjd)
%% - Pobs : matrix containing L3 obs [L3_base L3_rover]
%% - PosSat : matrix containing satellite position [X_base,Y_base,Z_base,X_rover,Y_rover,Z_rover] (m)
%% - ElevSat : matrix of satellite elevation (rad) [elev_base elev_rover]
%% - Dtropo : matrix of tropospheric delay correction (m) [dtropo_base dtropo_rover]
%% - sat_const : vector of satellite constellation [const] (1 = GPS, 2 = Galileo, 3 = Glonass)
%% - sat_PRN : vector of satellite PRN
%% - amb_index : matrix containing ambiguity number corresponding to obs [amb_base amb_rover]
%% - amb0 : approx ambiguities
%% - X0_base : base station coordinates [X,Y,Z] (m)
%% - X0_rover : approx rover station coordinates [X,Y,Z] (m)
%% - pivot : name of pivot satellite
%% - interv_pivot : interval of pivot satellite validity [mjd_min mjd_max]
%%
%% Output : 
%% - X =  4201575.42207108                  : X (m) 
%% - Y =  189859.220123781                  : Y (m) -> position in WGS84
%% - Z =  4779064.66690559                  : Z (m)
%% - DDN =                                  : % lambda*DDN values followed by [Nkn Nin Nkj Nij]_index
%%  		-19.22386498230120   11  1 14 4
%% 			-15.95452240039919   12  2 14 4
%%  		 12.70798318845537   13  3 14 4
%% 		 	 -5.90341014094431   15  5 14 4
%%  		-19.57968631928084   16  6 14 4
%%  		-10.35766380618947   17  7 14 4
%% 			-31.17191101499186   18  8 14 4
%%  		-15.65983940772206   19  9 14 4
%%  		 21.10299452446713   20 10 14 4
%%
%% - sigma02 =  0.0802353773556112          : sigma^2 of compensation
%% - V =                                    : residuals
%%    
%%         0.0316077318090322
%%         0.3940406712322329
%%        -0.4492030152309354
%%        -0.6061949610298285
%%        -0.9160045866227051
%%         0.3671790137783546
%%         0.2741622094259810
%%         1.9515150612008085
%%         0.2818285022353706
%%        -0.2943740050123722
%%         0.0261052661627730
%%        -1.3678729988680374
%%         1.0696252070105885
%%         8.7037829506630118
%%        -0.0458391305469945
%%         5.9964534601530062
%%        -0.8890634697682742
%%        -0.7874324639867609
%%         0.3396343510511883
%%    
%% - Vnorm =                                : normalized residuals
%%
%%         0.0316077318090322
%%         0.3940406712322329
%%        -0.4492030152309354
%%        -0.6061949610298285
%%        -0.9160045866227051
%%         0.3671790137783546
%%         0.2741622094259810
%%         0.9515150612008085
%%         0.2818285022353706
%%        -0.2943740050123722
%%         0.0261052661627730
%%        -0.3678729988680374
%%         0.0696252070105885
%%         0.7037829506630118
%%        -0.0458391305469945
%%         0.9964534601530062
%%        -0.8890634697682742
%%        -0.7874324639867609
%%         0.3396343510511883
%%
%% - Qxx =                                  : Var_covar matrix
%%    
%%       Columns 1 through 3:
%%    
%%         7.784096396444699   0.199696638416061   5.020430321182305
%%         0.199696638416061   3.291506563493738   1.136079887133114
%%         5.020430321182306   1.136079887133115   9.376077643598707
%%         6.882761948234403   0.528876708679703   7.169268720781800
%%         0.277181628248869   0.679543030146928   1.434044601350368
%%    
%%       Columns 4 and 5:
%%    
%%         6.882761948234403   0.277181628248868
%%         0.528876708679703   0.679543030146928
%%         7.169268720781802   1.434044601350368
%%         8.328530219783001  -0.236338941603696
%%        -0.236338941603696   2.633270027178538
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sat_index = cell(size(sat_const,1),1);
for i = 1:size(sat_index,1)

	if (sat_const(i) == 1)
		const = 'G';
	elseif (sat_const(i) == 2)
		const = 'E';
	elseif (sat_const(i) == 3)
		const = 'R';
	end
	
	sat = sprintf('%s%02d',const,sat_PRN(i));
	
	sat_index{i,1} = sat;

end

[result] = calc_LS_phase(t,Pobs,PosSat,ElevSat,Dtropo,sat_index,amb_index,amb0,X0_base,X0_rover,pivots,interv_pivot);

X = result.X;
Y = result.Y;
Z = result.Z;

DDN = result.DDN;

V = result.V;
Vnorm = result.Vnorm;

sigma02 = result.sigma02;
Qxx = result.Qxx;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
