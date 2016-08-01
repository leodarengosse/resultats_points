function [stat]=calc_stat_indic(Qxx,X,Y,Z);
%% [stat]=calc_stat_indic(Qxx,X,Y,Z);
%% Confidence ellipsoid computation (1-sigma) from Qxx matrix estimated in calc_LS
%% Correlation matrix computation
%% DOP computation
%% If several dtr estimated, TDOP correspond to first epoch
%%
%% Beilin Jacques - DPTS - 2010-11-10
%% Clement Fontaine - 05-10-2013
%% 
%% Input :
%% - Qxx : variance matrix (4*4)
%% - X,Y,Z : estimated coordinates
%%
%% Output :
%% - stat : structure containing statistic indicators : 
%% 		- ell : structure containing error ellipsoid in the local frame
%%		 example :
%%  		 ell =
%%  		 {
%%   		  Az1 = 319.00     % Azimuth 1
%%   		  Ele1 = -80.00    % Elevation 1
%%    		 err1 = 24.71     % Semi-major axis 1
%%    		 Az2 = 54.00      % Azimuth 2
%%    		 Ele2 = -11.00    % Elevation 2
%%    		 err2 = 2.84      % Semi-major axis 2
%%    		 Az3 = 357.00     % Azimuth 3
%%    		 Ele3 = 17.00     % Elevation 3
%%   		  err3 = 3.35      % Semi-major axis 3
%%  		 }
%% 		- Qenu : variance matrix in local frame
%% 		- Corr_enu : correlation matrix in local frame
%% 		- GDOP,PDOP,HDOP,VDOP,TDOP
%%   
%%   if X, Y or Z == 0, no calculation and 0 or empty matrix returned
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (X ==0 || Y == 0 || Z == 0 || size(Qxx,1)<4)
	stat.ell.Az1 = 0;
	stat.ell.Ele1 = 0;     
	stat.ell.err1 = 0;
	stat.ell.Az2 = 0;
	stat.ell.Ele2 = 0;
	stat.ell.err2 = 0;
	stat.ell.Az3 = 0;
	stat.ell.Ele3 = 0;
	stat.ell.err3 = 0;
	stat.Qenu = [];
	stat.Corr_enu = [];
	stat.GDOP = 0;
	stat.PDOP = 0;
	stat.HDOP = 0;
	stat.VDOP = 0;
	stat.TDOP = 0;
	return
end

ell=[];
[nl,nc]=size(Qxx);

% geographic coordinates
[l,p,h] = tool_cartgeo_GRS80(X,Y,Z);

% topocentric matrix ECEF
M = diag(ones(nc,1));
M(1:3,1:3)=[ -sin(p)*cos(l) -sin(p)*sin(l) cos(p)   
    -sin(l)         cos(l)        0        
    cos(p)*cos(l)   cos(p)*sin(l) sin(p) ];

% variance propagation
Qenu=M*Qxx*M';
stat.Qenu = Qenu;

% standard deviation extraction
sE2 = Qenu(1,1);
sN2 = Qenu(2,2);
sU2 = Qenu(3,3);
st2 = Qxx(4,4);

% DOP
stat.GDOP = sqrt(sE2+sN2+sU2+st2);
stat.PDOP = sqrt(sE2+sN2+sU2);
stat.HDOP = sqrt(sE2+sN2);
stat.VDOP = sqrt(sU2);
stat.TDOP = sqrt(st2);


% eigen values and vectors computation
[V, LAMBDA]=eig(Qenu(1:3,1:3));

% Error ellipsoid computation
[stat.ell.Az1,stat.ell.Ele1]=eigvec2azele(V(:,1));
stat.ell.err1 = LAMBDA(1,1);
[stat.ell.Az2,stat.ell.Ele2]=eigvec2azele(V(:,2));
stat.ell.err2 = LAMBDA(2,2);
[stat.ell.Az3,stat.ell.Ele3]=eigvec2azele(V(:,3));
stat.ell.err3 = LAMBDA(3,3);



% Correltation matrix
Corr_enu = zeros(nl,nc);

for i=1:nc
    for j=i:nc
        if i==j
            Corr_enu(i,j) = 1;    
        else
            Corr_enu(i,j) = Qenu(i,j) / sqrt(Qenu(i,i)) / sqrt(Qenu(j,j));
            Corr_enu(j,i) = Corr_enu(i,j);
        end
    end
end

stat.Corr_enu = Corr_enu;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






function [Az,ele]=eigvec2azele(V)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    X = V(1);
    Y = V(2);
    Z = V(3);

    XY = sqrt(X^2+Y^2);

    Az = 2 * atan(X/(Y+XY));
    ele = atan(Z/XY);

    Az = mod (Az,2*pi);

    Az = round(Az * 200 / pi);
    ele = round(ele  * 200 / pi);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
