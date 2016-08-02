function [CPR,H] = construc_mat(NAV_header,NAV_data,RNX_header,RNX_data,Point,mjd,epoch)
    
    l = length(Point.List_PRN);
    Station = get_station(RNX_header);%%r�cup�ration des coordonn�es de la station
    for i=1:1:l
        %Calcul des coordonn�es des satellites
        Sat = get_sat_data(NAV_header,NAV_data,RNX_header,RNX_data,Point.List_PRN(i),mjd,epoch);
        if (isempty(Sat))

        else
            Sat = pos_sat(Sat,mjd); %% ajout des coordonn�es calcul�es au temps t
            P_Dist_Point = dist_cart(Sat,Point);
            R = dist_cart(Sat, Station); %calcul de la distance entre la station et le satellite
            CPR(i,1) = (R-Sat.P_Dist); %cr�ation du vecteur des corrections
            [H(i,1),H(i,2),H(i,3)] = calc_deriv(Point,Sat); %cr�ation de la matrices des observations
            
        end
    end
end