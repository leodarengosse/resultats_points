function Sat = get_sat_data(NAV_header,NAV_data,RNX_header,RNX_data,PRN,tmjd,epoch)
        
        Eph = get_ephemeris(NAV_header,NAV_data,'G',PRN,tmjd);
        c = 299792458;
    
        Obs = get_obs(RNX_header,RNX_data,'G',PRN,epoch);

        
        if (isempty(Eph) || isempty(Obs))
            str = sprintf('PRN %s introuvables dans les fichiers',num2str(PRN));
            Sat = [];
        else
            Sat = Eph;
            Sat.dte = corr_dte_nav(Eph,tmjd);
            Sat.dt_relat = corr_dtrelat_nav(Eph,tmjd);
            Sat.C1 =Obs.C1;
            Sat.P_Dist = Obs.C1 + c * Sat.dte + c * Sat.dt_relat;
            %Sat.P_Dist = Obs.C1
        end

end
