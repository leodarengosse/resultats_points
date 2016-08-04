format long g
more off
clear all

%%
%Initialisation des variables

fichier_points = 'Points_mesures/Points.txt';
fichier_points_connus ='Points_mesures/Points_connus.txt';

numero_connu = 'B';
numero_mesure = 'B1';

fichier_n = 'Rinex/igp0200k15.16n'; %10:00

%fichier_o = 'Rinex/igp0200k05.16o'; %A 10:05
fichier_o = 'Rinex/igp0200k08.16o'; %B 10:08


%fichier_n = 'Rinex/igp0200o10.16n'; %14h00

%fichier_o = 'Rinex/igp0200o08.16o'; %A 14:08
%fichier_o = 'Rinex/igp0200o10.16o'; %B 14:10


%% 
%Chargement des donn�es
Points_data = load_point(fichier_points); %%chargement du fichier de points mesurees avec la tablette
Points_connus = load_point(fichier_points_connus); %chargement du fichier de points connus

Point = get_data_point(Points_data,numero_mesure); %%matrice de donn�es, num�ro du point
Point_connu = get_data_point(Points_connus,numero_connu); 

%%
% Calcul des coordonnees g�ocentriques en coordonn�es
[Point.X,Point.Y,Point.Z] = tool_geocart_GRS80(Point.Long,Point.Lat,Point.Alt);%coordonn�es � traiter

%test des fonctions de conversions
[Point.Long1,Point.Lat1,Point.Alt1] = tool_cartgeo_GRS80(Point.X,Point.Y,Point.Z);

%Conversion du point connu en coordonn�es
[Point_connu.X,Point_connu.Y,Point_connu.Z] = tool_geocart_GRS80(Point_connu.Long,Point_connu.Lat,Point_connu.Alt);

[Result_init] = compare_points(Point, Point_connu)

%%Stockage des coordon�es approch�e donn�es par la tablette
[Point.X_mes,Point.Y_mes,Point.Z_mes] = tool_geocart_GRS80(Point.Long,Point.Lat,Point.Alt);

%%
% Chargement des donn�es d'�ph�m�rides et d'observations

[NAV_header,NAV_data]=load_rinex_n(fichier_n); %%fichier de navigation

%%temps utc pour les points donc r�cup�ration du leap seconds
lpsec = NAV_header.LEAP_SECONDS; 

[RNX_header,RNX_data]=load_rinex_o(fichier_o); %%fichier d'observations

%%
% Calcul du temps et de l'�poque 

t = Point.Time + lpsec; %%temps de la mesure utc en secondes convertit en temps gps
%t = Point.Time ; %%temps de la mesure utc en secondes convertit en temps gps

tpointgps = s1970_t(t); %cr�ation d'une structure gps time

str = sprintf('%d:%d:%d',tpointgps.hh, tpointgps.min,tpointgps.sec) %affichage de la date pour v�rifier

mjd = tpointgps.mjd;  %% r�cup�ration du jour julien
epoch = get_epoch_from_mjd(RNX_header,mjd);  %%r�cup�ration de l'�poque correspondante

%%
%Estimation des nouvelles coodonn�es par moindres carr�s

for i=1:1:5
    %Calcul des matrices pour les moindres carr�s
    [B,CPR,H] = construc_mat(NAV_header,NAV_data,RNX_header,RNX_data,Point,mjd,epoch);  
    B
    dX=inv(H'*H)*H'*B; %solution des moindres carr�s
    
    V = H*dX - B; %calcul du vecteur des r�sidus
    n = length(B);
    p = length(dX);
    sigma_2 = (V'*V)/(n-p);%facteur unitaire de variance
    Point.X = Point.X + dX(1);
    Point.Y = Point.Y + dX(2);
    Point.Z = Point.Z + dX(3);
    compare_points(Point, Point_connu); %%ca diverge
end

%%
% Conversion en coordonn�es g�ographiques
[Point.Long,Point.Lat,Point.Alt] = tool_cartgeo_GRS80(Point.X,Point.Y,Point.Z);

%%
% Affichage des r�sultats
Result_init;
Result_fin = compare_points(Point, Point_connu);
