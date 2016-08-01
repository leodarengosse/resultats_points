function [dr]=corr_dtropo_saast(P,T,H,h,zen)
%%  [dr]=corr_dtropo_saast(P,T,H,h,zen)
%%  calculates ZPD from Saastamoinen model
%%
%% Jacques Beilin - ENSG - 2014-06-11
%% Clement Fontaine 2013-10-31
%%
%% Input :
%% - P     : pressure (hPa)  
%% - T     : temperature (K)
%% - H     : humidity (%)
%% - h     : ellipsoid height (m)
%% - zen   : zenithal angle (rad)
%%
%% Output
%% - dr    : zenithal tropospheric delay (m)
%%
%% Reference : SUBROUTINE TROPOS(Z,HS,T,P,RH,MODEL,DR) from Bernese GPS Software
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    bcor = [ 1.156 1.006 0.874 0.757 0.654 0.563 ];

    if (h>4000) 
        h=4000;
    end

    if(h<0) 
        h=0;
    end
    
    p = P * (1 - 2.26e-5 * h)^5.225;
    t = T - h * 0.0065;
    rhum = H * exp( -6.396e-4 * h);
    rhum = min(rhum,100);

    % WATER VAPOR PRESSURE
    e = ( rhum / 100 ) * exp( -37.2465 + 0.213166 * t - 2.56908e-4 * t * t);

    % MODEL 1: SAASTAMOINEN

    % HEIGHT IN KM
    hl = h / 1000;
    if (hl < 0) 
        hl = 0;
    end
    if (hl > 4) 
        hl = 4;
    end
    i = floor(hl) + 1;

    % REFERENCE HEIGHT FOR LINEAR INTERPOLATION IN TABLE BCOR
    href = i - 1;
    b = bcor(i) + (bcor(i+1) - bcor(i))*(hl - href);
	
    dr = 0.002277 * (p + (1255 / t + 0.05) * e - b * (sin(zen) / cos(zen))^2) / cos(zen);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
