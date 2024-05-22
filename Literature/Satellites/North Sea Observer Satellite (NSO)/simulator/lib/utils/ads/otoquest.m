%**************************************************************************
% Optimal Two-Observation Quaternion Estimation Method (OTOQuEst)
%
% Written by 05gr833
%
% Parameters:
%      Input:
%            r1    = u(1:3)   Sun vector in ECI reference frame.
%            r2    = u(4:6)   Magnetic field vector in ECI reference frame.
%            b1    = u(7:9)   Sun vector in SBRF.  
%            b2    = u(10:12) Magnetic field vector in SBRF.
%            a1    = u(13)    Weight of the first vector set.
%            a2    = u(14)    Weight of the second vector set.
%            q_old = u(15:18) Previously determined quaternion. 
%
%     Output:
%            q_opt  The optimal quaternion describing the attitude.
%
%**************************************************************************
function q_opt = otoquest(u)

% Input
r1=[u(1);u(2);u(3)];
r2=[u(4);u(5);u(6)];
b1=[u(7);u(8);u(9)];
b2=[u(10);u(11);u(12)];
a1=u(13);
a2=u(14);
q_old=[u(15);u(16);u(17);u(18)];

% Handling collinearity
if (cross(r1,r2)==0 | cross(b1,b2)==0)
    q_opt=q_old;
else

    % Third vector set is calculated, these are orthonormal vectors
    r3=cross(r1,r2)/norm(cross(r1,r2));
    b3=cross(b1,b2)/norm(cross(b1,b2));

    % Calculation of constants, which are used for determining the optimal
    % quaternion
    alpha=(1 + dot(b3,r3))*(a1*dot(b1,r1) + a2*dot(b2,r2))+dot(cross(b3,r3),(a1*cross(b1,r1) + a2*cross(b2,r2)));
    beta=dot((b3+r3),(a1*cross(b1,r1)+a2*cross(b2,r2)));
    gamma=sqrt(alpha^2+beta^2);

    % Computing the optimal attitude quaternion
    if alpha >= 0
        tmp  = [(gamma+alpha)*cross(b3,r3)+beta*(b3+r3);(gamma+alpha)*(1+dot(b3,r3))];
        q_opt=(1/(2*sqrt(gamma*(gamma+alpha)*(1+dot(b3,r3)))))*tmp;
    else
        tmp  = [beta*cross(b3,r3)+(gamma-alpha)*(b3+r3);beta*(1+dot(b3,r3))];
        q_opt=(1/(2*sqrt(gamma*(gamma-alpha)*(1+dot(b3,r3)))))*tmp;
    end

    % Handling quaternion 'flips'
    err=qdiff(q_old,q_opt);
    if err >= pi/2;
        q_opt=-q_opt;
    end

end

%**************************************************************************
% The qdiff function is used to find the difference between two quaternion
% (error), which is used for correcting the direction. q_err=q_att*q_ref^-1
function angle = qdiff(q1, q2)
q2(1:4) = [ -q2(1:3);
    q2(4) ];
A= [ q1(4) q1(3) -q1(2) q1(1);
    -q1(3) q1(4) q1(1) q1(2);
    q1(2) -q1(1) q1(4) q1(3);
    -q1(1) -q1(2) -q1(3) q1(4) ];
q = A*q2;
angle = 2*acos(q(4));