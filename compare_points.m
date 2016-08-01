function [Result] = compare_points(PointA, PointB)
    
    %Calcul des variations en chaque axe
    dx = PointA.X - PointB.X;%meters
    dy = PointA.Y - PointB.Y;%meters
    dz = PointA.Z - PointB.Z;%meters
    
%     scatter3(PointA.X,PointA.Y,PointA.Z,'red')
%     scatter3(PointB.X,PointB.Y,PointB.Z,'blue')
    
    dist = dist_cart(PointA,PointB); 
    
    %Calcul des variations
    dl = PointA.Lat - PointB.Lat; %degrees
    dL = PointA.Long - PointB.Long; %degres
    dh = PointA.Alt - PointB.Alt; %meters
    
    %Création d'une structure présentant les résultats
    Result = struct('dx', dx, 'dy', dy, 'dz', dz, 'Dist', dist, ...
                    'dl', dl, 'dL', dL, 'dh', dh);