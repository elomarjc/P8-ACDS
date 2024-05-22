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


% Persistent variables:
    % Two variables are declared to hold values for the next iteration.
    persistent x_hat_old    % old state estimate
    persistent P_old        % covarinace matrix
    
    % the rest is in order to speed things up... these are initialized once
    persistent bias_size
    persistent No_elements_Q
    persistent No_elements_z
    persistent No_elements_x_hat
    persistent sensor_model
    persistent Q
    persistent R
    persistent R_mean
    persistent W
    persistent I
    persistent x_tilde
    persistent Use_eksperimental_Q_and_R_matrices
    persistent Use_stored_state_info
    persistent Use_old_description
    persistent store_p_step

% Function inputs:
    B_Rsun_meas = input(1:3);   % Sun vector measurement in SBRF
    B_b_meas    = input(4:6);   % Magnetometer vector measurement in SBRF
    omega_meas  = input(7:12);  % Gyro vector measurement in SBRF - Format: [x+ y+ z+ x- y- z-]
   
    I_Rsun      = input(13:15); % Sun vector from onboard model in ECI
    I_b         = input(16:18); % magnetometer vector from onboard model in ECI
      
    h_mw_meas   = input(19:21); % Momentum wheel measurement in SBRF
    N_mt_meas   = input(22:24); % Magnetorquer measurement in SBRF
   
    sim_time    = input(25);    % Simulation time, - used to determine when to sample q0 and start the filter
    q0          = input(26:29); % Initial state of q
   
    Ts          = input(30);    % Sample time of the filter 
   
    Eclipse     = input(31);    % boolean input: 1 = eclipse, 0 = not eclipse
   

% The program:
  if sim_time < 2         % to change the startup delay change here... sim_time is in sec.
  % Initialization: 
      %   - This code will run for the first 2 seconds of the simulation. 
      %     It just feeds q0 and w_meas directly to the output.
        
        warning off % disable warnings in workspace
      
        store_p_step = 0;
        
      % Data selection:
        Use_stored_state_info = 1;  % 1 - using stored data to init P and bias states
                                    % 0 - used definitions in this file
        Use_old_description = 0;
                                    
                                    
        Use_eksperimental_Q_and_R_matrices = 1; % 1 - using experimental based data to init R and Q
                                                % 0 - used definitions in this file
        

        No_elements_Q = 9;  % size of Q referencing the system parameters, not incl. bias.
        No_elements_z = 12; % size of measurement vector
        bias_size     = 12;
        No_elements_x_hat = 10+bias_size;   % bias accounts for 6 states
        I = diag([1377 1623 1569]/1e6);     % Satellite inertia matrix
        
        if Use_eksperimental_Q_and_R_matrices
            load test/ads/ekf/kalman_R_matrix;
            R_mean = R;
        end
                                    
        if Use_stored_state_info == 1
            % bias startup state:
            
            load test/ads/ekf/kalman_bias_estimat.mat
            bias_startup_state = bias;
            %bias_startup_state(4:6) = [0 0 0];
            %bias_startup_state(1:3) = bias_startup_state(1:3)*5;
            
            % The old Bias_start_state
            if Use_old_description
                bias_startup_state = [0 0 0 0 0 0 0.00075980068378 0.00114684473875 0.00102070453291 -0.00076130017879 -0.00114420895840 -0.00099636058800];
            end
            
            % P startup state:
            P_old(10:15,10:15) = 1e-8*eye(6,6);
            P_old = zeros(21);
            
            %**
            load long_with_albedo/P_old_600.mat            
            P_old(1:3,1:3) = 1e-1*eye(3);
            %**
            
            if Use_old_description
                load stored_P_old
                P_old = stored_P_old;
            end
            
        else
            % bias init values = 0;
            bias_startup_state = zeros(1, bias_size);
            
            % P-matrix init
            P_a = [1e-3*eye(3) zeros(3,6);zeros(3) 1e-3*eye(3) zeros(3);zeros(3,6) 1e-3*eye(3)]; 
            P_b = zeros(No_elements_Q,bias_size);
            P_c = zeros(bias_size,No_elements_Q);
            P_d = [zeros(6,bias_size); zeros(6,6) 1e-8*eye(6,6)];
            %P_d = [2e-5*eye(3,3) zeros(3,9);zeros(3,3) 0*eye(3,3) zeros(3,6); zeros(6,6) 1e-8*eye(6,6)];
            P_old = [P_a P_b ; P_c P_d]; 
        end
        
       
        % Q-matrix setup - is not changed after this...
           %    q                               w                               h_mw
        if Use_eksperimental_Q_and_R_matrices
            load test/ads/ekf/kalman_Q_matrix
            Q_a = Q;
            Q_b = zeros(No_elements_Q,bias_size);
            Q_c = zeros(bias_size,No_elements_Q);

            % for simulation 1,2:
            Q_d = [1e-1*eye(6,bias_size); zeros(6,6) 1e-20*eye(6,6)]; 

            % for simulation 4:
            Q_d = [1e-1*eye(3,bias_size); zeros(3,3) 1e-3*eye(3,bias_size-3); zeros(6,6) 1e-20*eye(6,6)];
        else
            Q_a = [1e-4*eye(3) zeros(3,6);zeros(3) 1e-4*eye(3) zeros(3);zeros(3,6) 1e-4*eye(3)]; % Noise model for the predictor
            Q_b = zeros(No_elements_Q,bias_size);
            Q_c = zeros(bias_size,No_elements_Q);
            Q_d = [zeros(6,bias_size); zeros(6,6) 1e-8*eye(6,6)]; 
        end

        Q   = [Q_a Q_b ; Q_c Q_d];

        
        % init the state vector
        omega_meas  = (omega_meas(1:3) - omega_meas(4:6))/2;  % the two sets of omega measurements have opposite sign....
        x_hat_old   = [q0' omega_meas' h_mw_meas' bias_startup_state]';

        output = [x_hat_old(1:4); x_hat_old(5:7); zeros(12,1); bias_startup_state']; % just pass through the input
  else
  % This code will run after the startup delay

        % Test for eclipse. 
        % If in eclipse the sun-sensor-measurements are considered
        % unreliable, - the R-weigths are increased
        
        if Use_eksperimental_Q_and_R_matrices
            if (Eclipse == 1)
                %R(1:3,1:3) = 1e10*eye(3);
            else
                R = R_mean;
            end
        else
            if (Eclipse == 1) %| (sim_time > 2000)
                   % Sun Vector: 1e10                Mag. vector                       gyro                            gyro
                R = [1e10*eye(3) zeros(3,9);zeros(3) 1e-7*eye(3) zeros(3,6);zeros(3,6) 1e-6*eye(3) zeros(3);zeros(3,9) 1e-6*eye(3)]; % Noise model for the measurement
            else    %1e1
                R = [1e1*eye(3) zeros(3,9);zeros(3) 1e-7*eye(3) zeros(3,6);zeros(3,6) 1e-6*eye(3) zeros(3);zeros(3,9) 1e-6*eye(3)]; % Noise model for the measuremen
            end
        end

    % Get measurements and previous states:
    
        % The measurement vector z consists of
        % z = input(1:No_elements_z) = [B_Rsun_meas B_b_meas omega_meas]
        
        % Previous state vector and covariance matrix
        % x_hat_old;
        % P_old; 

    % Prediction: 
            % 1. Compute Phi, if R and Q are variable the put them here 
            Phi = jacobianF(x_hat_old,I,Ts);              % Jacobian of system matrix A

            % 2. Propagate the covariance matrix:
            P   = Phi*P_old*Phi' + Q;   % it is presumed that the noise i Gaussian, zero mean => W = eye(3)
                                        % else:  P = Phi*P_old*Phi' + WQW';

            % 3. Calculate the next estimate:
            step  = Ts;  % sample time
            Nrung = 50;  % number of steps during sample time
            Ndist = 0;   % disturbances are not modelled onboard
            x_hat = runge(Ndist, N_mt_meas, h_mw_meas, x_hat_old, I, step, Nrung);

            % 4 a). Compute estmated measurements 
            B_b_hat     = qmul( qmul(qinv(x_hat(1:4)), [I_b' 0]), x_hat(1:4)' );    % rotate I_b_meas from I frame to B frame
            B_b_hat     = B_b_hat(1:3)';
            B_Rsun_hat  = qmul( qmul(qinv(x_hat(1:4)), [I_Rsun' 0]), x_hat(1:4)' ); % rotate I_s_meas from I frame to B frame
            B_Rsun_hat  = B_Rsun_hat(1:3)';

            % 4 b).  Compute residual (Estimation Error)
            % Use this only when no quaternions in z:
            B_Rsun_tilde = B_Rsun_meas     - B_Rsun_hat - x_hat(11:13);
            B_b_tilde    = B_b_meas        - B_b_hat    - x_hat(14:16);        
            w_tilde_a    = omega_meas(1:3) - x_hat(5:7) - x_hat(17:19);
            w_tilde_b    = omega_meas(4:6) + x_hat(5:7) - x_hat(20:22);
            residual     = [B_Rsun_tilde ; B_b_tilde ; w_tilde_a; w_tilde_b];

   % Correction:
        % Compute Kalman Gain:

        % Create H matrix:
        H_a = [eye(3); -eye(3)];
        H_b = [2*skew_matrix(B_Rsun_hat)    zeros(3)    zeros(3)    ;
               2*skew_matrix(B_b_hat)       zeros(3)    zeros(3)    ; 
               zeros(6,3)                   H_a         zeros(6,3) ];

        H   = [H_b eye(12)];

        % Compute Kalman Gain
        K = P*H'*inv(H*P*H' + R);

        % Update estimate
        x_hat_new = addstates(x_hat, expandstate(K*(residual),x_hat(4))); 

        % Update Covariance Matrix
        P_new = (eye(No_elements_x_hat-1) - K*H)*P;       

        % Store internal states
        P_old = P_new;
        x_hat_old = x_hat_new;

        % Output result as a vector
        output = [x_hat_new(1:4); x_hat_new(5:7); residual(1:12); x_hat_old(11:22)];
  
  
%         if store_p_step < sim_time
%             store_p_step = store_p_step + 100;
%             save(sprintf('P_old_%d',store_p_step/100),'P_old');
%         end
  end

%**************************************************************************
%
%  The following section contains the functions used by the Kalman filter above
%

function F = jacobianF(x,I,Ts)

size_x = length(x-1);
size_bias = size_x - 10;

w    = x(5:7);
h_mw = x(8:10);
bias = x(11:size_x);

S_w     = skew_matrix(w);
Iw      = I*w;
S_Iw    = skew_matrix(I*w);
S_h_mw  = skew_matrix(h_mw);

Fsys_a  = [ -S_w                 .5*eye(3)                     zeros(3)   ;
             zeros(3)            inv(I)*(S_Iw-S_w*I + S_h_mw) -inv(I)*S_w ;
             zeros(3)            zeros(3)                      zeros(3)   ];       

Fsys    = [Fsys_a zeros(9,size_bias); zeros(size_bias,9) zeros(size_bias,size_bias)]; % add bias terms      
            
Gsys    = [zeros(3)     zeros(3);
 	       inv(I)       -inv(I) ;
           zeros(3)     eye(3)  ;
           zeros(12,3)  zeros(12,3)];

ss_sys  = ss(Fsys,Gsys,eye(length(Fsys)),0);    % Specify state-space model
dss_sys = c2d(ss_sys,Ts);                       % Ts - Sample time in sec.
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
function x = expandstate(x0, q4)
%expands the quaternion state from 3 elements to 4
%q4 of the old state q4
%q = x0(1:3);
%w = x0(4:6);

%expanding the quaternion states

qnorm = sqrt(1+(x0(1)^2) + (x0(2)^2) + (x0(3)^2));
qexp = [x0(1)/qnorm; x0(2)/qnorm; x0(3)/qnorm; 1/qnorm];

%Returns the expanded states
x = [qexp; x0(4:length(x0))];

%***************************************************************************
%function y = norm_q(x)
%q = x0(1:4);

%qnorm = sqrt((x(1)^2) + (x(2)^2) + (x(3)^2) + (x(4)^2));
%qexp = [x(1)/qnorm; x(2)/qnorm; x(3)/qnorm; x(4)/qnorm];
%norm(qexp);

%Return normalized state
%y = [qexp; x(5:length(x))];

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
    y  = y + (1/6)*(k1+2*k2+2*k3+k4);
end

qnorm  = norm(y(1:4));
y(1:4) = y(1:4)./qnorm;
result = y;

function rungres = f(x, I, Next, Nmw)
q        = x(1:4);
w        = x(5:7);
hmw      = x(8:10);
bias     = x(11:length(x));

q_dot    = 0.5.*[-skew_matrix(w) w; -w' 0]*q;
w_dot    = inv(I)*(Next + Nmw - skew_matrix(w)*(I*w + hmw));
hmw_dot  = Nmw;
bias_dot = zeros(length(bias),1);

rungres  = [q_dot; w_dot; hmw_dot; bias_dot];

%**************************************************************************
function xsub = subtractstates(x0);

size = length(x0)/2;

qx = x0(1:4);
xrest = x0(5:size);
qz    = x0(1+size:4+size);
zrest = x0(5+size:size+size);

qx = [-qx(1:3) x0(4)];
qnew = qmulkalman(qz,qx);

xtemp(1:4) = qnew;
xtemp(5:size) = zrest - xrest;
xsub = xtemp(1:size)';


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

size    = length(qhatk);

qnew    = qmulkalman(qK(1:4),qhatk(1:4)');

xtemp(1:4)    = qnew(1:4);
xtemp(5:size) = qhatk(5:size) + qK(5:size);
x             = [xtemp(1:4)'; xtemp(5:size)'];

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