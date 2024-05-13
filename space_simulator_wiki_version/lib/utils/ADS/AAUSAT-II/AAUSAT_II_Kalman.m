%**************************************************************************
% Kalman Filter for the AAUSAT-II using quaternion as input
%
% Written by group 833
%
% Parameters:
%   Input:  
%     q_meas      = input(1:4)  - Quaternion
%     omega_meas  = input(5:7)  - Gyro vector measurement in SBRF 
%                                    - Format: [x y z]
%
%     h_mw_meas   = input(19:21) - Momentum wheel measurement in SBRF
%     N_mt_meas   = input(22:24) - Magnetorquer measurement in SBRF
% 
%     sim_time    = input(25)    - Simulation time, - used to determine
%                                   when to sample q0 and start the filter
%     q0          = input(26:29) - Initial state of q
% 
%     Ts          = input(30)    - Sample time of the filter 
% 
%     Eclipse     = input(31)    - boolean input: 1 = eclipse,
%                                                 0 = not eclipse
% 
%   Output:
%     States:
%       output(1:4)     = q;      Filtered attitude
%       output(5:7)     = omega;  Filtered angular velocity
%
%     Residuals: (Kalman filter residual: measurement - estimate)
%       output(8:11)    = quaternion_tilde   
%       output(12:14)   = omega_tilde;
%
%**************************************************************************
function output = AAUSAT_II_kalman(input)

% Persisten variables
  % - variables are declared to hold values to the next iteration...
  persistent x_hat_old    % old state estimate
  persistent P_old        % covarinace matrix
  
  % the rest is in order to speed things up... these are initialized once
  persistent No_elements_z
  persistent No_elements_x_hat
  persistent sensor_model
  persistent Q
  persistent R
  persistent W
  persistent Estimat
  persistent I

% Function input:
  q_meas      = input(1:4);      %quaternion measurement
  omega_meas  = input(5:7);      %angular velocity measurement
  h_mw_meas   = input(8:10);     %Angular momentum measurement
  N_mt_meas   = input(11:13);    %Torque measurement
   
  sim_time    = input(14);
  q0          = input(15:18);
   
  Ts          = input(19);
  eclipse     = input(20);
   
% Initialization
  if sim_time < 5
    Ts = 0.1;   % Simulation stepsize
    No_elements_z = 10;
    No_elements_x_hat = 10;
    sensor_model = eye(No_elements_x_hat);
    % System data:
    I = diag([1377 1623 1569]/1e6);

    Q   = [1e-7*eye(3) zeros(3,6);zeros(3) 1e-6*eye(3) zeros(3);zeros(3,6) 1e-4*eye(3)]; % Noise model for the predictor

    x_hat_old   = [q0' omega_meas' h_mw_meas']'; 
    P_old       = 0.00001 * eye(No_elements_x_hat-1);
    output      = [x_hat_old(1:7);zeros(9,1)];
  else
    % test for eclipse:
    if eclipse ==1
      %If in eclipse do not trust quaternion
              % q                               w                    h_mw
      R   = [1e10*eye(3) zeros(3,6); zeros(3) 1e-9*eye(3) zeros(3);zeros(3,6) 1e-9*eye(3)]; % Noise model for the measurements
    else
      R   = [1e-3*eye(3) zeros(3,6); zeros(3) 1e-9*eye(3) zeros(3);zeros(3,6) 1e-9*eye(3)]; % Noise model for the measurements
    end
    
    z = input(1:No_elements_z);   % measurements

    % Predictors: 
      % 1. Compute Phi, if R and Q are variable the put them here 
      Phi = jacobianF(x_hat_old,I,Ts);              % Jacobian of system matrix A
                      
      % 2. Propagate the covariance matrix:
      % OK
      P   = Phi*P_old*Phi' + Q;   % it is persumed that the noise i Gausian, zero mean => W = eye(3)
                                  % else:  P = Phi*P_old*Phi' + WQW';
      % 3. Calculate the next estimat:
      step = Ts;  % sample time
      Nrung = 5;  % number of steps during sample time
      Ndist = 0;  % disturbances are not modelled on board
      x_hat = runge(Ndist, N_mt_meas,h_mw_meas , x_hat_old, I, step, Nrung);

      % 4 a). Compute estmated measurements vector z
      z_hat = (sensor_model*x_hat)';
      H = eye(No_elements_x_hat-1);

      % 4 b).  Compute residual (Estimation Error)
      % residual = z_hat' - z;
      residual = subtractstates([z_hat' z]);

    % Correction:
      % 5. Compute Kalman Gain:
      K = P*H'*inv(H*P*H'+ R);

      % 6. Update estimate
      x_hat_new = addstates(x_hat, expandstate(K*(remove_q4(residual))));

      % 7. Update Covariance Matrix
      P_new = (eye(No_elements_x_hat-1) - K*H)*P;       

      % store internal states
      P_old = P_new;
      x_hat_old = x_hat_new;

      % output result as a vector
      output = [x_hat_new(1:7);remove_q4(residual)];
  end
  
%**************************************************************************
%
%  The following section contains the functions used by the Kalman filter above
%
function F = jacobianF(x,I,Ts)

w    = x(5:7);
h_mw = x(8:10);

S_w     = skew_matrix(w);
Iw      = I*w;
S_Iw     = skew_matrix(I*w);
S_h_mw   = skew_matrix(h_mw);

Fsys    = [ -S_w         .5*eye(3)                     zeros(3)  ;
            zeros(3)    inv(I)*(S_Iw-S_w*I + S_h_mw) -inv(I)*S_w;
            zeros(3)    zeros(3)                      zeros(3) ];

Gsys    = [zeros(3) zeros(3);
 	       inv(I)   -inv(I) ;
           zeros(3) eye(3) ];



ss_sys  = ss(Fsys,Gsys,eye(9),0);           % Specify state-space model
dss_sys = c2d(ss_sys,Ts);                   % Ts - Sample time in sec.
PHIsys  = dss_sys.A;
F       = [PHIsys];

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
function x = expandstate(x0)
%expands the quaternion state from 3 elements to 4

%expanding the quaternion states

qnorm = sqrt(1+(x0(1)^2) + (x0(2)^2) + (x0(3)^2));
qexp = [x0(1)/qnorm; x0(2)/qnorm; x0(3)/qnorm; 1/qnorm];

%Returns the expanded states
x = [qexp; x0(4:9)];

%***************************************************************************

function result = runge(Ndist, Nmt,Nmw, x0, I, step, Nrung)
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
    y = y+(1/6)*(k1+2*k2+2*k3+k4);
end

qnorm = norm(y(1:4));
y(1:4)=y(1:4)./qnorm;
result = y;

function rungres = f(x, I, Next, Nmw)
q = x(1:4);
w = x(5:7);
hmw = x(8:10);
q_dot = 0.5.*[-skew_matrix(w) w; -w' 0]*q;
w_dot = inv(I)*(Next - Nmw - cross(w, ((I*w) + hmw)));
hmw_dot = Nmw;


rungres = [q_dot; w_dot; hmw_dot];

%**************************************************************************
function xsub = subtractstates(x0);

qx = x0(1:4);
xrest = x0(5:10);
qz = x0(11:14);
zrest = x0(15:20);

qx = [-qx(1:3) x0(4)];
qnew = qmulkalman(qz,qx);

xtemp(1:4) = qnew;
xtemp(5:10) = zrest(1:6) - xrest(1:6);
x2(1:10) = xtemp(1:10);
xsub = x2';

%**************************************************************************
function q = qmulkalman(q1, q2)

A=  [q1(4) q1(3) -q1(2) q1(1);
    -q1(3) q1(4) q1(1) q1(2);
     q1(2) -q1(1) q1(4) q1(3);
    -q1(1) -q1(2) -q1(3) q1(4)];

q = (A*q2')';

%**************************************************************************
function y = remove_q4(x)

y = [x(1:3); x(5:10)];

%**************************************************************************
function x = addstates(qhatk, qK) %(x0)

x0     = [qhatk qK];
qhatk  = x0(1:4);
qK     = x0(11:14);

qnew    = qmulkalman(qK,qhatk);

xtemp(1:4)  = qnew(1:4);
xtemp(5:10) = x0(5:10) + x0(15:20);
x           = [xtemp(1:4)'; xtemp(5:10)'];
