function [B,CPR,H,PosSat,PR,X0] = construc_mat(NAV_header,NAV_data,RNX_header,RNX_data,Point,mjd,epoch)
    
    l = length(Point.List_PRN);
    Station = get_station(RNX_header);%%r�cup�ration des coordonn�es de la station
    X0 = zeros(3,1);
    X0(1,1) = Station.X; 
    X0(2,1) = Station.Y;
    X0(3,1) = Station.Z;
    for i=1:1:l
        %Calcul des coordonn�es des satellites
        Sat = get_sat_data(NAV_header,NAV_data,RNX_header,RNX_data,Point.List_PRN(i),mjd,epoch);
        if (isempty(Sat))

        else
            Sat = pos_sat(Sat,mjd); %% ajout des coordonn�es calcul�es au temps t
            PosSat(i,1) = Sat.X;
            PosSat(i,2) = Sat.Y;
            PosSat(i,3) = Sat.Z;
            PR(i,1) = Sat.P_Dist;
          
            %distance entre les coordonn�es approch�es fix�es et les satellites
            Dist_Point_Sat_mes = dist_cart_coord(Point.X_mes,Point.Y_mes,Point.Z_mes,Sat.X,Sat.Y,Sat.Z); 
            
            %calcul de la distance entre la station et le satellite
            Dist_Station_Sat = dist_cart(Sat, Station); 
                         
% 
% S.dr = D[i] + epoch.cdtr -S.PR
% ou D : distance g�om�trique
% cdtr = biais d'horloge recepteur converti en distance
% PR : pseudo distance corrig�e


            %calcul de la distance entre le point et le satellite
            Dist_Point_Sat = dist_cart(Point,Sat); 
            
            %cr�ation du vecteur des corrections (diff�rence entre la
            %distance calcul�e entre le satellite et la station et la
            %pseudo-distance mesur�e
            
            CPR(i,1) = (Dist_Station_Sat + 125.717574997686 -Sat.P_Dist); 
            
            %Construction de la matrice B PRi(fixe) - CPR - R0(change � chaque it�ration)
            B(i,1) = Dist_Point_Sat_mes + CPR(i,1) - Dist_Point_Sat;
            
            %cr�ation de la matrice des observations
            [H(i,1),H(i,2),H(i,3)] = calc_deriv(Point,Sat); 
            
        end
    end
end