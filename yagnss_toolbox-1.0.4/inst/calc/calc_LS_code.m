function [result] = calc_LS_code(t,PosSat,Dobs,ElevSat,sat_index,X0,options)
%% function [result] = calc_LS_code(PosSat,Dobs,ElevSat,sat_index,X0,options) 
%% Least Square computation for code
%% GPS | GLONASS | GALILEO | GPS + GLONASS + GALILEO | GPS + GLONASS | GPS + GALILEO
%% Estimation of 1 cdtr per epoch, cGGTO and cGPGL for all data in input -> if one epoch
%% Estimation of one pos per epoch or one pos for the whole obs -> if several epochs
%%
%% GGTO = GPS to Galileo Time Offset
%% GPGL = GPS to GLonass time offset
%%
%%
%% P depends on elevation of satellites
%% 
%% Clement Fontaine 2013-12-18
%%
%%
%% Input : 
%% - t : vector containing time of obs (one time per obs) [t] (mjd)
%% - PosSat : matrix containing satellite position [X,Y,Z] (m)
%% - Dobs : vector containing observations (m)
%% - ElevSat : vector of satellite elevation (rad)
%% - sat_index : cell with satellite id {'constPRN'} . Format A1I2: ex {'G12';'G14';'G02'}
%% - X0 : initial values [X,Y,Z,cdtr,cGGTO,cGPGL]
%% - options : structure (optional)
%%     {
%%			constraint : postition estimation constraint in meters. If not defined or equal to 0: no constraint
%%     }
%%
%% Output : 
%% - result : structure containing results
%%    result =
%%    {
%%      t = 56442.0                            : mjd
%%      X =  4201575.42207108                  : X (m) 
%%      Y =  189859.220123781                  : Y (m) -> position in WGS84
%%      Z =  4779064.66690559                  : Z (m)
%%      cdtr =  108.729331571262               : cdtr (m)
%%      cGGTO = 0                              : c * GPS to Galileo Time Offset (m)
%%      cGPGL = -215.541766073757              : c * Glonass to GPSTime (m)
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
%%      nb_GPS =  10                           : number of GPS satellites used in compensation
%%      nb_GLO =  9                            : number of Glonass satellites used in compensation
%%      nb_GAL = 0                             : number of Galileo satellites used in compensation
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
%% If one dtr / epoch and one global pos : cdtr array in output, corresponding
%% to t in output 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%~ 
%~ PosSat(:,1).^2 + PosSat(:,2).^2 + PosSat(:,3).^2
%~ input('r')


%%%%% Set options

% Default parameters
constraint = 0;

% Options
if (nargin == 7)

	if isfield(options,'constraint') 
		if options.constraint~=0
			constraint = options.constraint;
		else
			constraint = 0;
		end
	else
		options.constraint = constraint;
	end
	
end

%%%%% Computation parameters : 

nb_GAL_GGTO = 2; % number of min Galileo satellites to compute GGTO 
% -> if less satellites than nb_GAL_GGTO, no Galileo satellites used in computation
nb_GLO_GPGL = 2; % number of min Glonass satellites to compute GPGL
% -> if less satellites than nb_GLO_GPGL, no Glonass satellites used in computation

%%%%% number of epochs
t_epoch = unique(t);  % unique returns elements without repetition
N_epoch = length(t_epoch);

%%%%% Approx values : 
% X0 = [X; Y; Z; cdtr; cGGTO; cGPGL];
if length(X0)<6
	X0 = [X0(:);zeros(6-length(X0),1)];
end

%%%%% Output -> if no computation, return input
result.t = 0;
result.X = X0(1);
result.Y = X0(2);
result.Z = X0(3);
result.cdtr = X0(4);
result.cGGTO = X0(5);
result.cGPGL = X0(6);
result.sigma02 = 0;
result.V = [];
result.Vnorm = [];
result.n_iter = 0;
result.Qxx = [];
result.nb_GPS = 0;
result.nb_GLO = 0;
result.nb_GAL = 0;
result.sat_index = [];

sta_pos = X0(1:3); % keep initial coord if position is constraint
pos0 = X0(1:3); 

cdtr = X0(4);
cGGTO = X0(5);
cGPGL = X0(6);

X0 = [];



%%%%% Test number of satellites

% GPS satellites
nb_GPS = sum(cell2mat(strfind(sat_index(:,1),'G')));
% Glonass satellites
nb_GLO = sum(cell2mat(strfind(sat_index(:,1),'R')));
% Galileo satellites
nb_GAL = sum(cell2mat(strfind(sat_index(:,1),'E')));
% Total satellite number
nb_sat = nb_GPS + nb_GLO + nb_GAL;

% No computation if just Gal + Glo
if (nb_GPS == 0 && nb_GLO > 0 && nb_GAL > 0)
	tool_print_info('No computation with only GAL and GLO satellites',3);
	return
end


%%%%% number of satellites to compute a position

% only one constellation -> no estimation of any time offset
if (nb_GPS == nb_sat || nb_GLO == nb_sat || nb_GAL == nb_sat)

	estim_GPGL = 0;
	estim_GGTO = 0;
	
	% number of satellites to compute a position
	nb_sat_min = 4; % min 4 satellites to compute (position + dtr)
	
	X0 = [X0;pos0;cdtr*ones(N_epoch,1)];
	pos_cdtr = size(X0,1)-N_epoch+1;
	
% several constellations -> time offsets estimation
else 

	% number of satellites to compute a position
	nb_sat_min = 6;  % 3 for position + 1 for dtr + 1 for GGTO + 1 for GPGL = 6
	
	X0 = [X0;pos0];
	
	% cGGTO
	if (nb_GAL<nb_GAL_GGTO) % too few Galileo satellites to estimate the offset between GPSTime and GalileoTime
	
		estim_GGTO = 0;
		nb_GAL = 0;
		nb_sat_min = nb_sat_min - 1;
		GAL_pos = cell2mat(sat_index(:,1))=='E';
		
		% suppress Galileo satellites
		PosSat(GAL_pos==1,:)=[];
		Dobs(GAL_pos==1)=[];
		ElevSat(GAL_pos==1)=[];
		sat_index(GAL_pos==1,:) = [];
		t(GAL_pos==1,:) = [];
		
	else
		estim_GGTO = 1;
		X0 = [X0;cGGTO];	
		pos_cGGTO = size(X0,1);
		
	end
	
	% cGPGL
	if (nb_GLO<nb_GLO_GPGL) % too few Glonass satellites to estimate the offset between GPSTime and GlonassTime

		estim_GPGL = 0;
		nb_GLO = 0;
		nb_sat_min = nb_sat_min - 1;
		GLO_pos = cell2mat(sat_index(:,1))=='R';
		
		% suppress GLONASS satellites
		PosSat(GLO_pos==1,:)=[];
		Dobs(GLO_pos==1)=[];
		ElevSat(GLO_pos==1)=[];
		sat_index(GLO_pos==1,:) = [];
		t(GLO_pos==1,:)=[];
		
	else
	
		estim_GPGL = 1;
		X0 = [X0;cGPGL];
		pos_cGPGL = size(X0,1);
		
	end
	
	% cdtr
	X0 = [X0;cdtr*ones(N_epoch,1)];
	pos_cdtr = size(X0,1)-N_epoch+1;

	
end

%%%%% Total satellite number
nb_sat = nb_GPS + nb_GLO + nb_GAL;

% nb of required satellites to inverse matrix
if nb_sat <= nb_sat_min
	if nb_GAL>=4
		tool_print_info(sprintf('Only GAL , nb_GAL == 4'),2);% debug
	else
		tool_print_info(sprintf('%2d satellites < %2d required satellites : no computation',nb_sat,nb_sat_min),1);
	return
	end
end

n = nb_sat; % equation number
p = length(X0); % unknown number

%%%%% Satellite positions
Xs=PosSat(:,1);
Ys=PosSat(:,2);
Zs=PosSat(:,3);



%%%%% Weight matrix
sigma = 2;
%~ P = eye(n)./4;
P = diag( (cos(pi/2-ElevSat)./sigma).^2 );
%~ P = eye(n);
% for glonass P = P/2 (we consider Var(glo) = Var(gps) / sqrt(2)
P_glo = find(cell2mat(sat_index(:,1))=='R');
P(P_glo,P_glo) = P(P_glo,P_glo)/2;

%%%%%% position constraints
if (constraint ~= 0)
	P2 = zeros(n+3,n+3);
	P2(1:n,1:n) = P;
	P2(n+1,n+1) = 1/(constraint^2);
	P2(n+2,n+2) = 1/(constraint^2);
	P2(n+3,n+3) = 1/(constraint^2);
	P = P2;
end

%%%%% position of GLO and GAL sat
pos_Glonass = (cell2mat(sat_index(:,1))=='R');
pos_Glonass = pos_Glonass(:,1);
pos_Galileo = (cell2mat(sat_index(:,1))=='E');
pos_Galileo = pos_Galileo(:,1);


%%%%% Convergence parameters and variable initialization
n_iter = 0;
sigma02priori = 1E25;
epsilon = 1E-6;

sigma02 = 0;
V = [];
Qxx = [];

while (n_iter<15) 

	% approx value of D
	D = ( (X0(1) - Xs).^2 + (X0(2) - Ys).^2 + (X0(3) - Zs).^2).^0.5;

	% estimated or constraint position
	if constraint == 0
		A = zeros(n,p);
		A(:,1:3) = [ (X0(1) - Xs)./D (X0(2) - Ys)./D (X0(3) - Zs)./D ];
		B = Dobs - D;
	else
		A = zeros(n+3,p);
		A(1:n,1:3) = [ (X0(1) - Xs)./D (X0(2) - Ys)./D (X0(3) - Zs)./D ];
		A(n+1:n+3,1:3) = diag(ones(3,1));
		B = [Dobs-D;sta_pos-X0(1:3)];
	end
        
    % cdtr
	for nb_dtr_estim = 1:N_epoch
		pos_dtr_estim = (t == t_epoch(nb_dtr_estim));
		A(1:n,pos_cdtr+nb_dtr_estim-1) = pos_dtr_estim;
		B(1:n) = B(1:n) - X0(pos_cdtr+nb_dtr_estim-1)*pos_dtr_estim;
    end
    
    % cGPGL
    if estim_GPGL==1 
		A(1:n,pos_cGPGL) = pos_Glonass;
		B(1:n) = B(1:n) - X0(pos_cGPGL)*pos_Glonass;
	else
		B(1:n) = B(1:n);
	end
    
    % cGGTO  
    if estim_GGTO==1 
		A(1:n,pos_cGGTO) = pos_Galileo;
		B(1:n) = B(1:n) - X0(pos_cGGTO)*pos_Galileo;
	else
		B(1:n) = B(1:n);
	end

    % system resolution
 
    %~ Qxx = inv(A' * P * A);
	% inversion after cholesky decomposition (conditionning pb with inv)
    U = chol(A'*P*A);
    invU = inv(U);
    Qxx = invU*invU';
        
    D = (A' * P * B);

    dX = Qxx * D;
    X0 = X0 + dX;
        
    % residual computation
    V = B - A*dX;
       
	if(n-p~=0) % cas GAL==4 debug
		sigma02 = (V' * P * V)/ (n-p);
	else
		sigma02 = 1000;
	end

	
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

if(n-p~=0) % cas GAL==4 debug
	sigma02 = (V' * P * V)/ (n-p);
else
	sigma02 = 1000;
end

Vnorm = V.*sqrt(diag(P));
Qxx = Qxx*sigma02;

% if constraint ~= 0, supp V corresponding to constraint

if constraint ~= 0
	V = V(1:end-3);
	Vnorm = Vnorm(1:end-3);
end

% output

result.X = X0(1);
result.Y = X0(2);
result.Z = X0(3);
result.cdtr = X0(pos_cdtr:pos_cdtr+N_epoch-1);

if estim_GGTO == 1
	result.cGGTO = X0(pos_cGGTO);
else
	result.cGGTO = NaN;
end

if estim_GPGL == 1
	result.cGPGL = X0(pos_cGPGL);
else
	result.cGPGL = NaN;
end 

result.t = unique(t);
result.sat_index = sat_index;
result.sigma02 = sigma02;
result.V = V;
result.Vnorm = Vnorm;
result.n_iter = n_iter;
result.Qxx = Qxx;

result.nb_GPS = sum(cell2mat(strfind(sat_index(:,1),'G')));
result.nb_GLO = sum(cell2mat(strfind(sat_index(:,1),'R')));
result.nb_GAL = sum(cell2mat(strfind(sat_index(:,1),'E')));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
