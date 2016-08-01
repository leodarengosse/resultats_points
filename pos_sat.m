function [Sat] = pos_sat(Sat,mjd)

    a = (Sat.sqrt_a)^2; %demi-grand axe
   
    GM = 3.986005e14; %geocentric gravtitational constant
    
    we = 7.2921151467e-5; %earth rotation
    
    n0 = sqrt(GM/(a^3));%mean motion

	% 1er pb : differentes echelles de temps
    tgps = mjd_t(mjd); % initialisation d'une structure gpstime à parir du mjd
    tk = tgps.wsec - Sat.TOE; %%GPS system time

    n = n0 + Sat.delta_n; %corrected mean motion
    
    Mk = Sat.M0 + n*tk; %mean motion corrected

    e = Sat.e;%eccentricity

    Ek = kepler_E(e,Mk);%eccentric anomaly

    sin_nu_k = (sqrt(1-e^2)*sin(Ek))/(1-e*cos(Ek));%true anomaly
    
    % 2 ere erreur  : il fautr le sin et le cos pour calculer nu_k car cela varie en 0 et 2*pi
    cos_nu_k = (cos(Ek) - e) / (1-e*cos(Ek));
    
    nu_k = atan2(sin_nu_k,cos_nu_k); % on fait atan2 pour regler le pb des 4 quadrants
    
    % nu_k = asin(sin_nu_k)
    
    %%nu_k = 2*atan(sqrt((1+e)/(1-e))*tan(Ek/2))
    
    phi_k = nu_k + Sat.omega; %argument of latitude
    
    delta_uk = Sat.cuc*cos(2*phi_k) + Sat.cus*sin(2*phi_k); %argument of latitude corrections

    delta_rk = Sat.crc*cos(2*phi_k) + Sat.crs*sin(2*phi_k); %radius corrections

    delta_ik = Sat.cic*cos(2*phi_k) + Sat.cis*sin(2*phi_k); %Inclination corrections

    uk = phi_k + delta_uk; %corrected argument of latitude

    rk = a*(1-e*cos(Ek)) + delta_rk; %Corrected radius

    ik = Sat.i0 + Sat.IDOT*tk + delta_ik; %corrected inclination
    
    %Position in the orbital Plane    
    Xk_prim = rk*cos(uk) ;
    Yk_prim = rk*sin(uk);
    
    %Corrected longitude of ascending node
    OMEGA_k = Sat.OMEGA + (Sat.OMEGA_DOT - we)*tk - we * Sat.TOE;
    
   % OMEGA_k = Sat.OMEGA + (Sat.OMEGA_DOT)*tk
    
    %Earth fixed geocentric satellite coordinates

    Sat.X = Xk_prim*cos(OMEGA_k) - Yk_prim*sin(OMEGA_k)*cos(ik);
    Sat.Y = Xk_prim*sin(OMEGA_k) + Yk_prim*cos(OMEGA_k)*cos(ik);
    Sat.Z = Yk_prim*sin(ik);
    
    % 3eme pb : il semble manquer la rotation de la terre
%     OMEGAe = -7.2921151467e-5 * tgps.wsec;    
%     [Sat.X,Sat.Y,Sat.Z] = tool_rotZ(Sat.X,Sat.Y,Sat.Z,OMEGAe);
%     
  
