function [dr_dx, dr_dy, dr_dz] = calc_deriv(Point,Sat)

    dr_dx = (Point.X-Sat.X)/dist_cart(Point,Sat);
    dr_dy = (Point.Y-Sat.Y)/dist_cart(Point,Sat);
    dr_dz = (Point.Z-Sat.Z)/dist_cart(Point,Sat);
    
    