%**************************************************************************
% Kalman Filter for the AAUSAT-II
% 
% Written by group 833
%
% Parameters:
%   Input:  
%     B_Rsun_meas = input(1:3)   - Sun vector measurement in SBRF
%     B_b_meas    = input(4:6)   - Magnetometer vector measurement in SBRF
%     omega_meas  = input(7:12)  - Gyro vector measurement in SBRF 
%                                    - Format: [x+ y+ z+ x- y- z-]
% 
%     I_Rsun      = input(13:15) - Sun vector from onboard model in ECI
%     I_b         = input(16:18) - magnetometer vector from onboard model in ECI  
% 
%     h_mw_meas   = input(19:21) - Momentum wheel measurement in SBRF
%     N_mt_meas   = input(22:24) - Magnetorquer measurement in SBRF
% 
%     sim_time    = input(25)    - Simulation time, - used to determine when
%                                    to sample q0 and start the filter
%     q0          = input(26:29) - Initial state of q
% 
%     Ts          = input(30)    - Sample time of the filter 
% 
%     Eclipse     = input(31)    - boolean input: 1 = eclipse, 0 = not eclipse
% 
%   Output:
%     States:
%       output(1:4)     = q;      Filtered attitude
%       output(5:7)     = omega;  Filtered angular velocity
%     Residuals: (Kalman filter residual: measurement - estimate)
%       output(8:10)    = Rsun_tilde;   
%       output(11:13)   = b_tilde;
%       output(14:16)   = w_tilde+;
%       output(17:19)   = w_tilde-;
%
%**************************************************************************
function output = AAUSAT_II_kalman(input)
    
% Function inputs:
    x         = input(1:10); 
    h_mw_meas = input(8:10);
    N_mt_meas = input(11:13);
    Ndist     = input(14:16);
    Ts        = input(17);
    I         = diag(input(18:20));

    % 3. Calculate the next estimate:
    Nrung = 5;  % number of steps during sample time
    x_hat = runge(Ndist, N_mt_meas, h_mw_meas, x, I, Ts, Nrung);

    output = x_hat;
%**************************************************************************
%
%  The following section contains the functions used by the Kalman filter above
%

%***************************************************************************
function output = skew_matrix(x)
% Returns a skew symmetric matrix of the input vector.

if length(x) ~= 3
    error('skew_matix: dimension error, - input is not a 1x3 vector')
end
    
output =  [0    -x(3)  x(2);
           x(3)  0    -x(1);
          -x(2)  x(1)  0  ];
                        
%***************************************************************************

function result = runge(Ndist, Nmt, Nmw, x0, I, step, Nrung)
% Step is the time to jump forward
% Nrung is the number of runs to do it in

Next = Ndist + Nmt;

y=x0;
stepsize = step/Nrung;

for kn=1:Nrung,
    k1 = stepsize*f(y       , I, Next, Nmw);
    k2 = stepsize*f(y+0.5*k1, I, Next, Nmw);
    k3 = stepsize*f(y+0.5*k2, I, Next, Nmw);
    k4 = stepsize*f(y+k3    , I, Next, Nmw);
    y  = y + (1/6)*(k1+2*k2+2*k3+k4);
end

qnorm  = norm(y(1:4));
y(1:4) = y(1:4)./qnorm;
result = y;

function rungres = f(x, I, Next, Nmw)
q        = x(1:4);
w        = x(5:7);
hmw      = x(8:10);
%bias     = x(11:length(x));

q_dot    = 0.5.*[-skew_matrix(w) w; -w' 0]*q;
w_dot    = inv(I)*(Next + Nmw - skew_matrix(w)*(I*w + hmw));
hmw_dot  = Nmw;
%bias_dot = zeros(length(bias),1);

rungres  = [q_dot; w_dot; hmw_dot];%; bias_dot];

%**************************************************************************
function q = qmulkalman(q1, q2)

A=  [q1(4) q1(3) -q1(2) q1(1);
    -q1(3) q1(4) q1(1) q1(2);
     q1(2) -q1(1) q1(4) q1(3);
    -q1(1) -q1(2) -q1(3) q1(4)];

q = (A*q2')';

%**************************************************************************
function y = qinv(q)

y = [-q(1:3); q(4)];

%**************************************************************************
function qres = qmul(q1, q2)

qres1 =  q1(1)*q2(4) + q1(2)*q2(3) - q1(3)*q2(2) + q1(4)*q2(1);
qres2 = -q1(1)*q2(3) + q1(2)*q2(4) + q1(3)*q2(1) + q1(4)*q2(2);
qres3 =  q1(1)*q2(2) - q1(2)*q2(1) + q1(3)*q2(4) + q1(4)*q2(3);
qres4 =  q1(4)*q2(4) - q1(1)*q2(1) - q1(2)*q2(2) - q1(3)*q2(3);

qres = [qres1 qres2 qres3 qres4];