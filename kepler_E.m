function Ek = kepler_E(e, M)
% œ
% This function uses Newton's method to solve Kepler's
% equation E - e*sin(E) = M for the eccentric anomaly,
% given the eccentricity and the mean anomaly.
%
% E - eccentric anomaly (radians)
% e - eccentricity, passed from the calling program
% M - mean anomaly (radians), passed from the calling program
% pi - 3.1415926...
%
% User M-functions required: none
% ------------------------------------------------------------

    %...Set an error tolerance:
    error = 1.e-8;
    %...Select a starting value for E:
    
    Ek=M;
    
    %...Iterate on Equation 3.14 until E is determined to within
    %...the error tolerance:
    delta_Ek = 1;
    while abs(delta_Ek) > error
        delta_Ek = (M + e*sin(Ek) - Ek)/(1 - e*cos(Ek));
        Ek=Ek+delta_Ek;
    end
    % 

end

%


