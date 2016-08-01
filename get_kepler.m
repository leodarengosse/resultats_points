function Sat = get_sat_data(NAV_header,NAV_data,PRN,tmjd)
 
Sat.Pdist 
    for i=1:1:length(List_PRN)
        
        if (isempty(Eph))
            sprintf('PRN %s introuvable',num2str(List_PRN(i)))
        else
            Sat = get_ephemeris(NAV_header,NAV_data,'G',List_PRN(i),tmjd)

            [Obs] = get_obs(RNX_header,RNX_data,'G',PRN,epoch);
            P_dist = 
        end
    end