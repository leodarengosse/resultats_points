function [B,CPR,H] = construc_mat(NAV_header,NAV_data,RNX_header,RNX_data,Point,mjd,epoch)
    
    l = length(Point.List_PRN);
    Station = get_station(RNX_header);%%récupération des coordonnées de la station
    for i=1:1:l
        %Calcul des coordonnées des satellites
        Sat = get_sat_data(NAV_header,NAV_data,RNX_header,RNX_data,Point.List_PRN(i),mjd,epoch);
        if (isempty(Sat))

        else
            Sat = pos_sat(Sat,mjd); %% ajout des coordonnées calculées au temps t
            
            %distance entre les coordonnées approchées fixées et les satellites
            Dist_Point_Sat_mes = dist_cart_coord(Point.X_mes,Point.Y_mes,Point.Z_mes,Sat.X,Sat.Y,Sat.Z); 
            
            %calcul de la distance entre la station et le satellite
            Dist_Station_Sat = dist_cart(Sat, Station); 
            
            %calcul de la distance entre le point et le satellite
            Dist_Point_Sat = dist_cart(Point,Sat) 
            
            %création du vecteur des corrections (différence entre la
            %distance calculée entre le satellite et la station et la
            %pseudo-distance mesurée
            
            CPR(i,1) = (Dist_Station_Sat-Sat.P_Dist); 
            
            %Construction de la matrice B PRi(fixe) - CPR - R0(change à chaque itération)
            B(i,1) = Dist_Point_Sat_mes + CPR(i,1) - Dist_Point_Sat;
            
            %création de la matrice des observations
            [H(i,1),H(i,2),H(i,3)] = calc_deriv(Point,Sat); 
            
        end
    end
end