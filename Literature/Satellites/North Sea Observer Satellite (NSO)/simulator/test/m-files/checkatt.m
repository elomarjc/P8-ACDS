function [OK q_ref] = checkatt(q_i,q_r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function checks whether the calculated reference quaternion rotates
% the angular momentum bias vector outside the first octant. If it does the
% value false is returned otherwise it returns true.
%
%   Input
%           q_i     Initial attitude quaternion
%           q_r     Reference attiude quaternion
%
%   Output
%           OK      Value is either 1 or 0
%           q_ref   Reference attitude quaternion for which the rotated
%                   angular momentum bias vector does not leave the first 
%                   octant.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OK=true;
q_ref = q_r;

h_bias_att = [1 1 1]'/sqrt(3);
h_bias_ext = [h_bias_att; 0];

q_i_conj = [-q_i(1:3) q_i(4)];
q_err = qmulatt(q_r,q_i_conj);
q_err_conj = [-q_err(1:3) q_err(4)];

h_bias_rot = qmulatt(q_err,h_bias_ext');
h_bias_rot = qmulatt(h_bias_rot,q_err_conj);

for p=1:3
    if(sign(h_bias_rot(p)) < 0)
        OK=false;
        
        theta = rand*(pi/4);
        phi   = rand*(pi/4);
        alpha = rand*(pi/4);
        
        q_x =[sin(theta/2) 0 0 cos(theta/2)];
        q_y =[0 sin(phi/2) 0 cos(phi/2)];
        q_z = [0 0 sin(alpha/2) cos(alpha/2)];
                
        q_ref = qmulatt(q_i,q_x);
        q_ref = qmulatt(q_ref,q_y);
        q_ref = qmulatt(q_ref,q_z);
    end
end


