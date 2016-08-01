function [result] = calc_LS_phase(t,Pobs,PosSat,ElevSat,Dtropo,sat_index,amb_index,amb0,X0_base,X0_rover,pivot,interv_pivot)
%% function [result] = calc_LS_phase(t,Pobs,PosSat,ElevSat,Dtropo,sat_index,amb_index,amb0,X0_base,X0_rover,pivot,interv_pivot)
%% Rover position and ambiguity estimation for phase processing
%% P depends on elevation of satellites
%%
%% Clement Fontaine 2013-12-18
%%
%% Input : 
%% - t : vector containing time of obs (one time per obs) [t] (mjd)
%% - Pobs : matrix containing L3 obs [L3_base L3_rover]
%% - PosSat : matrix containing satellite position [X_base,Y_base,Z_base,X_rover,Y_rover,Z_rover] (m)
%% - ElevSat : matrix of satellite elevation (rad) [elev_base elev_rover]
%% - Dtropo : matrix of tropospheric delay correction (m) [dtropo_base dtropo_rover]
%% - sat_index : cell with satellite id {'constPRN'} . Format A1I2: ex {'G12';'G14';'G02'}
%% - amb_index : matrix containing ambiguity number corresponding to obs [amb_base amb_rover]
%% - amb0 : approx ambiguities
%% - X0_base : base station coordinates [X,Y,Z] (m)
%% - X0_rover : approx rover station coordinates [X,Y,Z] (m)
%% - pivot : name of pivot satellite
%% - interv_pivot : interval of pivot satellite validity [mjd_min mjd_max]
%%
%% Output : 
%% - result : structure containing results
%%    result =
%%    {
%%      t =                                    : t (unique(t))
%%			56442.0
%%          56442.1
%%          56442.2
%%          56442.3

%%      X =  4201575.42207108                  : X (m) 
%%      Y =  189859.220123781                  : Y (m) -> position in WGS84
%%      Z =  4779064.66690559                  : Z (m)
%%      DDN =                                  : % lambda*DDN values followed by [Nkn Nin Nkj Nij]_index
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
%%      sigma02 =  0.0802353773556112          : sigma^2 of compensation
%%      V =                                    : residuals
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
%%      Vnorm =                                : normalized residuals
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
%%      n_iter =  2                            : number of iterations
%%      Qxx =                                  : Var_covar matrix
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
%%      nb_GPS =  150                           : number of GPS satellites obs used in compensation
%%      nb_GLO =  90                            : number of Glonass satellites obs used in compensation
%%      nb_GAL = 0                              : number of Galileo satellites obs used in compensation
%%      sat_index =
%%      
%%        {
%%          [1,1] = G01
%%          [2,1] = G03
%%          [3,1] = G06
%%          [4,1] = G11
%%          [5,1] = G14
%%          [6,1] = G19
%%          [7,1] = G20
%%          [8,1] = G22
%%          [9,1] = G28
%%          [10,1] = G32
%%          [11,1] = R06
%%          [12,1] = R07
%%          [13,1] = R08
%%          [14,1] = R13
%%          [15,1] = R14
%%          [16,1] = R21
%%          [17,1] = R22
%%          [18,1] = R23
%%          [19,1] = R32  
%%        }
%%      
%%    }
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% constants 
c = 299792458.0;

%~ f1 = 1575.42;
%~ f2 = 1227.60;
%~ lambda1 = c/(f1*1e6);
%~ lambda = lambda1 * (f1^2/(f1^2-f2^2));

sigma = 0.02; % sigma_l3

%output
result.t = 0;
result.X = X0_rover(1);
result.Y = X0_rover(2);
result.Z = X0_rover(3);
result.sigma02 = 0;
result.V = [];
result.Vnorm = [];
result.n_iter = 0;
result.Qxx = [];
result.DDN = [];
result.nb_GPS = 0;
result.nb_GLO = 0;
result.nb_GAL = 0;
result.sat_index = cell(0);


if(size(Pobs,1)<7)
	tool_print_info('No enough obs to compute a position : min = 7',3);
	return
end

% For PosSat -> use orbits computed with base data (for dtr : pos fixed)
% Matrix initialization
% for each epoch : nb_sat visible sat = nb_sat - 1 double differences < size(Pobs,1)
DDPobs = zeros(size(Pobs,1),1);
DDtropo = zeros(size(Pobs,1),1);
DDRconst = zeros(size(Pobs,1),1);
DDN_index = []; % each line = corresponding DDN
DDN = []; % DDN values followed by [Nkn Nin Nkj Nij]
Sigma_DD = zeros(size(Pobs,1),1);
DDSat = cell(size(Pobs,1),2);


% mjd of each epoch
t_ep = unique(t);

nb_obs = 1;

% Approximated values and double differences computation
% Set DDN_index + number of DDN to estimate

for ep = 1:size(t_ep,1)

	% position of local obs
	index_ep = find(t == t_ep(ep));
		
	% local obs 
	Pobs_ep = Pobs(index_ep,:);
	PosSat_ep = PosSat(index_ep,:);
	ElevSat_ep = ElevSat(index_ep,:);
	Dtropo_ep = Dtropo(index_ep,:);
	sat_index_ep = sat_index(index_ep,:);
	amb_index_ep = amb_index(index_ep,:);
		
	% is pivot visible ? 
	% current pivot
	current_pivot = find((interv_pivot(:,1)<=t_ep(ep)) .*(interv_pivot(:,2)>t_ep(ep)));

	pos_pivot = find(ismember(sat_index_ep(:,1),pivot{current_pivot}));
	
	if(pos_pivot<=0)
		% skip epoch
		continue;
	end

	% for each satellite
	for num_sat = 1:size(Pobs_ep,1)	
		
		if(num_sat ~= pos_pivot)
					
			% DDPobs
			% DDPobs = Pobs(rover,sat) - Pobs(base,sat) - Pobs(rover,pivot) + Pobs(base,pivot)
			DDPobs(nb_obs,1) = Pobs_ep(num_sat,2) - Pobs_ep(num_sat,1) - Pobs_ep(pos_pivot,2) + Pobs_ep(pos_pivot,1);

			% DDN index
			DDN0 = amb0(amb_index_ep(num_sat,2),2) - amb0(amb_index_ep(num_sat,1),1) - amb0(amb_index_ep(pos_pivot,2),2) + amb0(amb_index_ep(pos_pivot,1),1);
			DDN_loc = [amb_index_ep(num_sat,2) amb_index_ep(num_sat,1) amb_index_ep(pos_pivot,2) amb_index_ep(pos_pivot,1)];
			
			% new DDN ?
			if(size(DDN,1)==0)
				DDN = [DDN;DDN0,DDN_loc];
				DDN_index(nb_obs,1) = size(DDN,1);
			else
			
				pos_DDN = find(ismember(DDN(:,2:end),DDN_loc,'rows'));
				
				if length(pos_DDN)==0 % yes -> add new DDN in DDN_index
				
					DDN = [DDN;DDN0,DDN_loc];
					DDN_index(nb_obs,1) = size(DDN,1);
	
				else % no -> DDN already exists
				
					DDN_index(nb_obs,1) = pos_DDN;
					
				end
			end
			
			% DDtropo
			% DDtropo = DDtropo(rover,sat) - DDtropo(base,sat) - DDtropo(rover,pivot) + DDtropo(base,pivot)
			DDtropo(nb_obs,1) = Dtropo_ep(num_sat,2) - Dtropo_ep(num_sat,1) - Dtropo_ep(pos_pivot,2) + Dtropo_ep(pos_pivot,1);
				
			% DDR - constant part
			rin = sqrt( ( X0_base(1) - PosSat_ep(num_sat,1) )^2 + ( X0_base(2) - PosSat_ep(num_sat,2 ) )^2 + ( X0_base(3) - PosSat_ep(num_sat,3) )^2 ); % fixed
			rij = sqrt( ( X0_base(1) - PosSat_ep(pos_pivot,1) )^2 + ( X0_base(2) - PosSat_ep(pos_pivot,2) )^2 + ( X0_base(3) - PosSat_ep(pos_pivot,3) )^2 ); % fixed
			DDRconst(nb_obs,1) = - rin + rij; 
			
			% Sigma_DD
			% Sigma_DD^2 = sigma(rover,sat)^2 - sigma(base,sat)^2 - sigma(rover,pivot)^2 + sigma(base,pivot)^2
			Sigma_DD(nb_obs,1) = (sigma/cos(pi/2-ElevSat_ep(num_sat,2)))^2 + (sigma/cos(pi/2-ElevSat_ep(num_sat,1)))^2 + (sigma/cos(pi/2-ElevSat_ep(pos_pivot,2)))^2 + (sigma/cos(pi/2-ElevSat_ep(pos_pivot,1)))^2;
			
			
			% DDSat : index of current sat and pivot sat
			DDSat(nb_obs,1) = sat_index_ep(num_sat);
			DDSat(nb_obs,2) = sat_index_ep(pos_pivot);
			
			nb_obs = nb_obs +1;					
				
		end
		
	end
	
end
nb_obs = nb_obs - 1;

% redim matrix
DDPobs = DDPobs(1:nb_obs,:);
DDtropo = DDtropo(1:nb_obs,:);
DDRconst = DDRconst(1:nb_obs,:);
DDN_index = DDN_index(1:nb_obs,:);
Sigma_DD = Sigma_DD(1:nb_obs,:);

DDSat = DDSat(1:nb_obs,:);

%~ % debug
%~ for i = 1:size(DDN,1)
%~ 
%~ 
	%~ pos_sat = find(ismember(DDN_index,i));
	%~ 1234567
	%~ figure
	%~ DDN(i,:)
	%~ plot(DDPobs(pos_sat))
%~ 
%~ end
%~ % end debug


% Compensation

n_iter = 0;
sigma02priori = 1E25;
epsilon = 1E-6;



P= diag(1./Sigma_DD); % ???????????????
%~ P = diag(sigma^2*ones(nb_obs,1));


X0 = [X0_rover(1:3);DDN(:,1)];

while (n_iter<15) 

	A = zeros(nb_obs, 3 + size(DDN,1) ); % X Y Z + n*DDamb  
	B = zeros(nb_obs,1);
	
	
	pos = 1;
	for ep = 1:size(t_ep,1) % local obs
	
		% position of local obs
		index_ep = find(t == t_ep(ep));
			
		% local obs 
		Pobs_ep = Pobs(index_ep,:);
		PosSat_ep = PosSat(index_ep,:);
		amb_index_ep = amb_index(index_ep,:);
		sat_index_ep = sat_index(index_ep,:);

		% is pivot visible ? 
		% current pivot
		current_pivot = find((interv_pivot(:,1)<=t_ep(ep)) .*(interv_pivot(:,2)>t_ep(ep)));
		pos_pivot = find(ismember(sat_index_ep(:,1),pivot{current_pivot}));
	
		if(pos_pivot<=0)
			% skip epoch
			continue;
		end
		
		% for each satellite
		for num_sat = 1:size(Pobs_ep,1)	

			if(num_sat ~= pos_pivot)	
			
				% matrix A
					
				rkn = sqrt( ( X0(1) - PosSat_ep(num_sat,1) )^2 + ( X0(2) - PosSat_ep(num_sat,2) )^2  + ( X0(3) - PosSat_ep(num_sat,3) )^2 );
				rkj = sqrt( ( X0(1) - PosSat_ep(pos_pivot,1) )^2 + ( X0(2) - PosSat_ep(pos_pivot,2) )^2 + ( X0(3) - PosSat_ep(pos_pivot,3) )^2 );
				
				% dr/dx   dr/dy   dr/dz 
				A(pos,1) = ( X0(1) - PosSat_ep(num_sat,1) ) / rkn - ( X0(1) - PosSat_ep(pos_pivot,1) ) / rkj;
				A(pos,2) = ( X0(2) - PosSat_ep(num_sat,2) ) / rkn - ( X0(2) - PosSat_ep(pos_pivot,2) ) / rkj;
				A(pos,3) = ( X0(3) - PosSat_ep(num_sat,3) ) / rkn - ( X0(3) - PosSat_ep(pos_pivot,3) ) / rkj;
				
				% ambiguities
				A(pos, 3 + DDN_index(pos)) = - 1; %lambda;
								
				% matrix B
				DDR = rkn - rkj + DDRconst(pos);
				B(pos,1) = DDPobs(pos,1) - DDtropo(pos,1) - DDR + DDN(DDN_index(pos),1);
							
				pos = pos + 1;
				
			end
		
		end
	
	end
	
	% inversion after cholesky decomposition (conditionning pb with inv)
	
    U = chol(A' * P * A);
    invU = inv(U);
    Qxx = invU*invU';
        
    D = (A' * P * B);

    dX = Qxx * D;
    X0 = X0 + dX;
    
	DDN(:,1) = X0(4:end);
	        
    % residual computation
    V = B - A*dX;
       
	sigma02 = (V' * P * V)/ (nb_obs-(3+size(X0,1)));
	
	% convergence criteria (sigma0 variation < epsilon)
    n_iter = n_iter+1;
    if (abs((sigma02-sigma02priori))<epsilon)    		
        break
    end
    sigma02priori = sigma02; 
    
    % Normalized residuals
	Vnorm = V.*sqrt(diag(P));

	% change weight if Vnorm > 3
	err_pos = find(abs(Vnorm)>3);
	for err = 1:length(err_pos)
		P(err_pos(err),err_pos(err)) = P(err_pos(err),err_pos(err)) / abs(Vnorm(err_pos(err)));
	end
     
end

Vnorm = V.*sqrt(diag(P));
Qxx = Qxx*sigma02;

% results

result.t = t_ep;
result.X = X0(1);
result.Y = X0(2);
result.Z = X0(3);

result.sigma02 = sigma02;
result.V = V;
result.Vnorm = Vnorm;
result.n_iter = n_iter;
result.Qxx = Qxx;
result.DDN = DDN;
result.sat_index = sat_index;

% GPS satellites
result.nb_GPS = sum(cell2mat(strfind(sat_index(:,1),'G')));
% Glonass satellites
result.nb_GLO = sum(cell2mat(strfind(sat_index(:,1),'R')));
% Galileo satellites
result.nb_GAL = sum(cell2mat(strfind(sat_index(:,1),'E')));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
