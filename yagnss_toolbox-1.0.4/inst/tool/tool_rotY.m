function [Xr,Yr,Zr] = tool_rotY(X,Y,Z,alpha)
%% [Xr,Yr,Zr] = tool_rotY(X,Y,Z,alpha)
%%
%% Rotation around Y-axis
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
    
    R = [  cos(alpha)  0   sin(alpha)
           0           1   0
          -sin(alpha)  0   cos(alpha) ];
     
    B = R * [X ; Y ; Z];

    Xr = B(1);
    Yr = B(2);
    Zr = B(3);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
