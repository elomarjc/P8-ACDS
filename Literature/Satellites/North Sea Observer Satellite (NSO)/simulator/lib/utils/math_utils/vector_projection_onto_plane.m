%**************************************************************************
% vector_projection_onto_plane.m
%
% This function projects a vector v onto a plane with the given normal
% vector n.
%
% Written by 06gr1032
%
% Parameters:
%
%      Input:    
%            n = u(1:3) Normal vector to the plane.
%            v = u(4:6) Vector for projection. 
%
%     Output:
%            pro_vec    The projected vector on the plane.
%
%**************************************************************************
function pro_vec = vector_projection_onto_plane(u)

% Input
n = [u(1);u(2);u(3)];
v = [u(4);u(5);u(6)];

% Calculation of the projection matrix P.
P = (1/(n(1)^2 + n(2)^2 + n(3)^2))*[n(2)^2+n(3)^2 -n(1)*n(2) -n(1)*n(3);
                                -n(1)*n(2) n(1)^2+n(3)^2 -n(2)*n(3); 
                                -n(1)*n(3) -n(2)*n(3) n(1)^2+n(2)^2];

% Calculation of the projected vector.
pro_vec = P*v;
                            