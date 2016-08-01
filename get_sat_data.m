function Sat = get_sat_data(NAV_header,NAV_data,RNX_header,RNX_data,PRN,tmjd,epoch)
        
        Eph = get_ephemeris(NAV_header,NAV_data,'G',PRN,tmjd);
        Obs = get_obs(RNX_header,RNX_data,'G',PRN,epoch);
        
        if (isempty(Eph) || isempty(Obs))
            str = sprintf('PRN %s introuvables dans les fichiers',num2str(PRN));
            Sat = [];
        else
            Sat = Eph;
            Sat.P_Dist = Obs.C1;
        end

end
