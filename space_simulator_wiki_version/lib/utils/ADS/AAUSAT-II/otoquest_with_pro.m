%**************************************************************************
% Optimal Two-Observation Quaternion Estimation Method (OTOQuEst)
%
% Written by 05gr833
%
% Parameters:
%      Input:
%            r1       = u(1:3)   Sun vector in ECI.
%            r2       = u(4:6)   Magnetic field vector in ECI.
%            b1       = u(7:9)   Sun vector in SBRF.  
%            b2       = u(10:12) Magnetic field vector in SBRF.
%            a1       = u(13)    Weight of the first vector set.
%            a2       = u(14)    Weight of the second vector set.
%            w        = u(15:17) The angular velocity vector.
%            eclipse  = u(18)    If 1 then the satellite is in eclipse. 
%                                If 0 then the satellite is not in eclipse. 

%     Output:
%            q_opt  The optimal quaternion describing the attitude.
%
%**************************************************************************
function q_opt = otoquest(u)
persistent q_old;

% Input
r1=[u(1);u(2);u(3)];
r2=[u(4);u(5);u(6)];
b1=[u(7);u(8);u(9)];
b2=[u(10);u(11);u(12)];
a1=u(13);
a2=u(14);
w = u(15:17);
eclipse = u(18);

step  = 0.1;
Nrung=5;

if eclipse == 0
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
        if ~isempty(q_old)
            err=qdiff(q_old,q_opt);
            if err >= pi/2;
                q_opt=-q_opt;
            end
        else
            q_opt = -q_opt; 
            q_old = q_opt;
        end

        q_old = q_opt;
    end
else
    %during the first run q_old has to be initialized if Eclipse == 1...
    if isempty(q_old)
        q_old = [0 sqrt(0.91) 0 0.3]';
        disp('q_old is initialized...')
    end
    
    q_hat = runge(w, q_old, step, Nrung);

    q_old = q_hat;

    q_opt =q_hat;
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


%**************************************************************************
function result = runge(omega, q0, step, Nrung)
% Step is the time to jump forward
% Nrung is the number of runs to do it in

y=q0;
stepsize = step/Nrung;

for kn=1:Nrung,
    k1 = stepsize*f(y       ,omega);
    k2 = stepsize*f(y+0.5*k1,omega);
    k3 = stepsize*f(y+0.5*k2,omega);
    k4 = stepsize*f(y+k3    ,omega);
    y = y+(1/6)*(k1+2*k2+2*k3+k4);
end

qnorm = norm(y(1:4));
y(1:4)=y(1:4)./qnorm;
result = y;

function rungres = f(x,u)
q = x;
w = u;
q_dot = 0.5.*[-skew_matrix(w) w; -w' 0]*q;

rungres = q_dot;

%***************************************************************************
function output = skew_matrix(x)
% Returns a skew symmetric matrix of the input vector.

if length(x) ~= 3
    error('skew_matix: dimension error, - input is not a 1x3 vector')
end
    
output =  [0    -x(3)  x(2);
           x(3)  0    -x(1);
          -x(2)  x(1)  0  ];
