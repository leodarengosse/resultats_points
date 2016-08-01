function Dist = dist_cart(A,B)

    Dist = sqrt((A.X-B.X)^2+(A.Y-B.Y)^2+(A.Z-B.Z)^2);
end
