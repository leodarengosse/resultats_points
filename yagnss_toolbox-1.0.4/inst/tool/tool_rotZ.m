function [Xr,Yr,Zr] = tool_rotZ(X,Y,Z,alpha)
%% [Xr,Yr,Zr] = tool_rotZ(X,Y,Z,alpha)
%%
%% Rotation around Z-axis
%%
%% Beilin Jacques - ENSG/DPTS - 2012-05-16
%%
%% Input
%% - X,Y,Z : cartesian coordinates (m)
%% - alpha : rotation angle (rad)
%%
%% Output
%% - X,Y,Z : cartesian coordinates (m)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    
    R = [ cos(alpha)   -sin(alpha)   0
          sin(alpha)    cos(alpha)   0
          0             0            1];
     
    B = R * [X ; Y ; Z];

    Xr = B(1);
    Yr = B(2);
    Zr = B(3);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
