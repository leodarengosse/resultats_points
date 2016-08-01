function [X,Y,Z,dte,debug] = orb_from_eph(Eph,mjd)
%% function [X,Y,Z,dte,debug] = orb_from_eph(Eph,mjd)
%%
%% Orbit and satellite clock error calculation from ephemeris
%%
%% Clement Fontaine 2013-11-14
%%
%% Input :
%% - Eph : GPS or Galileo Eph structure, from get_ephemeris()
%%         Requested fields : 
%%                                - Eph.mjd
%%                                - Eph.TOE
%%                                - Eph.M0
%%                                - Eph.e
%%                                - Eph.delta_n
%%                                - Eph.sqrt_a
%%                                - Eph.omega
%%                                - Eph.cus
%%                                - Eph.cuc
%%                                - Eph.crs
%%                                - Eph.crc
%%                                - Eph.i0
%%                                - Eph.IDOT
%%                                - Eph.OMEGA
%%                                - Eph.OMEGA_DOT
%%                                - Eph.cic
%%                                - Eph.cis
%%                                - Eph.alpha0
%%                                - Eph.alpha1
%%                                - Eph.alpha2
%% - mjd : modified julian day
%%
%% Output : 
%% - X, Y, Z : satellite position in cartesian coordinates
%%              - WGS84 for GPS
%%              - GTRF for Galileo
%% - dte : satellite clock error
%% - debug : structure containing intermediate results
%%
%% Returns X = 0, Y = 0, Z = 0, dte = 0 and debug = cell if Eph.const 
%% different from 'G' or 'E'
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = 0;
Y = 0;
Z = 0;
dte = 0;
debug = cell(0);

if ~isfield(Eph,'const')
	return;
end

if (strcmp(Eph.const,'G') || strcmp(Eph.const,'E')) % GPS and Galileo

        % Earth orbit parameters (the same for GPS and Galileo)
        
        TOC = Eph.mjd; % Time of clock in mjd
        TOE = Eph.TOE; % Time of ephemeris in week seconds
        M0 = Eph.M0;
        e = Eph.e;
        Dn = Eph.delta_n;
        a = Eph.sqrt_a * Eph.sqrt_a;
        vieux_pi = Eph.omega; % vieux pi
        Cus = Eph.cus;
        Cuc = Eph.cuc;
        Crs = Eph.crs;
        Crc = Eph.crc;
        i0 = Eph.i0;
        iDOT = Eph.IDOT;
        OMEGA0 = Eph.OMEGA;
        OMEGA_DOT = Eph.OMEGA_DOT;
        Cic = Eph.cic;
        Cis = Eph.cis;

        tgps = mjd_t(mjd);
   
        if isfield(Eph,'gps_wk')
                gps_wk = Eph.gps_wk;
        elseif isfield(Eph,'gal_wk')
                gps_wk = Eph.gal_wk;
        else
                tool_print_info('gps_ws is not defined, position and clock error set to 0',2);
                return;
        end
        
        %~ if(abs(tgps.mjd-mjd)>) % If tgps is not close to mjd...
                %~ printf
                %~ return
        %~ end
        
        % Keplerian elements
        
        % Position on the orbit
        
        % Mean motion computation
        
        debug.mu = 398600440000000;
        debug.n0 = sqrt( debug.mu /a/a/a);
        debug.n = debug.n0 + Dn;
                
        % Mean anomaly computation at t
        
        t0 = (mjd - TOC) * 86400;
     
        debug.t0 = t0;
        M = M0 + debug.n * t0;
        debug.M = M;
        
        % Kepler's equation for eccentric anomaly E
        
        E0 = M;
        ct_debug=0;
        while 1==1
            E = M + e * sin(E0);
            if (abs(E-E0)<1e-12||ct_debug>20)
               break   
            end
            E0 = E;
            ct_debug=ct_debug+1;
        end
        debug.E = E;
        
        
        % True anomaly computation
        
        debug.v = 2 * atan((tan(E/2))*sqrt((1+e)/(1-e)));
        
        % Argument of latitude
        
        debug.phi = debug.v + vieux_pi;
        
        % Argument on latitude correction
        
        debug.dphi = Cus * sin(2*debug.phi) + Cuc * cos(2*debug.phi);
        
        % Radius correction
        
        debug.dr = Crs * sin(2*debug.phi) + Crc * cos(2*debug.phi);
        
        % Radius
        
        debug.r = a * (1 - e *cos(debug.E));

 
        % Position in orbital plane
        
        debug.xy = [ (debug.r + debug.dr) * cos(debug.phi + debug.dphi) 
                     (debug.r + debug.dr) * sin(debug.phi + debug.dphi)
                        0.0];
           
        % Position of orbital plane in the space         
      
        % temporal evolution af i and OMEGA  
        
        i = i0 + iDOT * debug.t0;
        OMEGA = OMEGA0 + OMEGA_DOT * debug.t0;
        debug.OMEGA = OMEGA;
        

        
        % inclination correction
        
        debug.di = Cic * cos(2*debug.phi) + Cis * sin(2*debug.phi); 
        debug.i = i + debug.di;
        
        % Keplerian elements to WGS84 cartesian coordinates
        
        
        [debug.ixy(1),debug.ixy(2),debug.ixy(3)] = tool_rotX(debug.xy(1),debug.xy(2),debug.xy(3),debug.i);
        
        [debug.X_ECI(1),debug.X_ECI(2),debug.X_ECI(3)] = tool_rotZ(debug.ixy(1),debug.ixy(2),debug.ixy(3),debug.OMEGA);
                
        % x,y in orbital plane
        % X in celestrian equatorial frame

        % Earth is spinning
        % omega_e * number of week seconds
        wsec = tgps.wsec;
        
        % if mjd<Eph.mjd and wk<Eph.gps_wk, wsec = - offset between mjd and Eph.mjd
        if (tgps.mjd<Eph.TOC && tgps.wk<gps_wk)
                wsec = wsec - 604800.0;
        end

        % if week of Eph < week of tgps, add 86400*7
        if (tgps.wk>gps_wk)        
                wsec = wsec + 604800.0;
        end
        
        
        OMEGAe = -7.2921151467e-5 * wsec;
        
        [debug.X_ECEF(1),debug.X_ECEF(2),debug.X_ECEF(3)] = tool_rotZ(debug.X_ECI(1),debug.X_ECI(2),debug.X_ECI(3),OMEGAe);
        
        % WGS84 or GTRF coordinates of SV antenna phase center position at time t
        
        X = debug.X_ECEF(1);
        Y = debug.X_ECEF(2);
        Z = debug.X_ECEF(3);   
           
                            
        % Relativist correction computation
        F = -4.442807633E-10;
        dt_relat = F * sqrt(a) * e * sin(E); 
        debug.dt_relat = dt_relat;
        
        % Satellite clock error computation
        dte =  Eph.alpha0 + Eph.alpha1 * t0 + Eph.alpha2 * t0^2; 
        debug.dte = dte;
       
        
else
        tool_print_info('Use orb_from_RK4 for Glonass position determination : coordinates and dte set to 0',2);
        return

end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
