function output = AAUSAT_II_Kalman_BS_vectors(input)

% To variables are declared to hold values to the next iteration...
persistent x_hat_old    % old state estimate
persistent P_old        % covarinace matrix

persistent No_elements_z
persistent No_elements_x_hat
persistent sensor_model
persistent Q
persistent R
persistent W
persistent x_tilde

% input to the function:
   B_Rsun_meas = input(1:3);
   B_b_meas    = input(4:6);   
   omega_meas  = input(7:9);
   I_Rsun      = input(10:12);
   I_b         = input(13:15);   
   
   h_mw_meas   = input(16:18);
   N_mt_meas   = input(19:21);
   
   sim_time    = input(22);
   q0          = input(23:26);
   
   Ts          = input(27);
   
   Eclipse     = input(28);
   
% Initialization: 
%   - This code will run for the first 0.5 seconds of the simulation. 
%     It just feeds q0 and w_meas directly to the output.

if sim_time < 0.5
    Ts = 0.1;       % Simulation stepsize
    No_elements_z = 9;
    No_elements_x_hat = 10;
    
    Q   = [1e-4*eye(3) zeros(3,6);zeros(3) 1e-9*eye(3) zeros(3);zeros(3,6) 1e-9*eye(3)]; % Noise model for the predictor
    R   = 0.0001 * eye(9);      % Noise model for the measurements
    %W   = eye(No_elements_z);  % omitted since W is estimated eye()

    x_hat_old   = [q0' omega_meas' h_mw_meas']'; %zeros(No_elements_x_hat,1);
    P_old       = 0.001 * ones(No_elements_x_hat-1);
    output      = [x_hat_old(1:7)];
else
    z = input(1:9);   % measurements    
 
    % Test for eclipse
    if Eclipse == 1
        Q = [1e4*eye(3) zeros(3,6);zeros(3) 1e-9*eye(3) zeros(3);zeros(3,6) 1e-9*eye(3)];
    else
        Q = [1e-4*eye(3) zeros(3,6);zeros(3) 1e-9*eye(3) zeros(3);zeros(3,6) 1e-9*eye(3)];
    end
    
    % System data:
    I = diag([1.3e-4 2e-4 1e-4]);

    % Previous state vector and covariance matrix
    % x_hat_old;
    % P_old; 

    % Predictors: 
        % 1. Compute Phi, if R and Q are variable the put them here 
        Phi = jacobianF(x_hat_old,I,Ts);              % Jacobian of system matrix A
                        
        % 2. Propagate the covariance matrix:
        P   = Phi*P_old*Phi' + Q;   % it is persumed that the noise i Gausian, zero mean => W = eye(3)
                                    % else:  P = Phi*P_old*Phi' + WQW';

        % 3. Calculate the next estimat:
        step  = Ts;  % sample time
        Nrung = 5;   % number of steps during sample time
        Ndist = 0;   % disturbances are not modelled on board
        x_hat = runge(Ndist, N_mt_meas, h_mw_meas, x_hat_old, I, step, Nrung);
        
        % 4 a). Compute estmated measurements vector z
        %C_q         = calculate_C_q(x_hat(1:4));
        
        B_b_hat     = qmul( qmul(qinv(x_hat(1:4)), [I_b' 0]), x_hat(1:4)' );    % rotate I_b_meas from I frame to B frame
        B_b_hat     = B_b_hat(1:3)';
        B_Rsun_hat  = qmul( qmul(qinv(x_hat(1:4)), [I_Rsun' 0]), x_hat(1:4)' ); % rotate I_s_meas from I frame to B frame
        B_Rsun_hat  = B_Rsun_hat(1:3)';
       
        % 4 b).  Compute residual (Estimation Error)
        % Use this with no quaternions in z:
        B_b_tilde    = B_b_meas - B_b_hat;
        B_Rsun_tilde = B_Rsun_meas - B_Rsun_hat;
        w_tilde      = omega_meas - x_hat(5:7);
        residual     = [B_b_tilde ; B_Rsun_tilde ; w_tilde];

    % Correction:
    % 5. Compute Kalman Gain:
    H = [2*skew_matrix(B_b_hat)     zeros(3)    zeros(3);
         2*skew_matrix(B_Rsun_hat)  zeros(3)    zeros(3); 
         zeros(3)                   eye(3)      zeros(3)]; 

    K = P*H'/(H*P*H' + R);

    % 6. Update estimate
    x_hat_new = addstates(x_hat, expandstate(K*(residual))); 

    % 7. Update Covariance Matrix
    P_new = (eye(No_elements_x_hat-1) - K*H)*P;       

    % store internal states
    P_old = P_new;
    x_hat_old = x_hat_new;

    % output result as a vector
    output = [x_hat_new(1:7)];
end
%***************************************************************************

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
%q = x0(1:3);
%w = x0(4:6);

%expanding the quaternion states

qnorm = sqrt(1+(x0(1)^2) + (x0(2)^2) + (x0(3)^2));
qexp = [x0(1)/qnorm; x0(2)/qnorm; x0(3)/qnorm; 1/qnorm];

%Returns the expanded states
x = [qexp; x0(4:9)];

%***************************************************************************
function y = norm_q(x)
%q = x0(1:4);

qnorm = sqrt((x(1)^2) + (x(2)^2) + (x(3)^2) + (x(4)^2));
qexp = [x(1)/qnorm; x(2)/qnorm; x(3)/qnorm; x(4)/qnorm];
norm(qexp);

%Return normalized state
y = [qexp; x(5:10)];

%***************************************************************************

function result = runge(Ndist, Nmt,Nmw, x0, I, step, Nrung)
% Step is the time to jump forward
% Nrung is the number of runs to do it in

Next = Ndist + Nmt;

y=x0;
%N_mw = x0(8:10); 
%Iscb = [1.3e-4 0 0; 0 2e-4 0; 0 0 1e-4]
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
hmw = Nmw;
%bias = x(11:13);
q_dot = 0.5.*[-skew_matrix(w) w; -w' 0]*q;
w_dot = inv(I)*(Next - Nmw - cross(w, ((I*w) + hmw)));
hmw_dot = Nmw;
%bias_dot = bias;


%rungres = [q_dot; w_dot; hmw_dot; bias_dot];
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

%**************************************************************************
function C_q = calculate_C_q(q)

C_q = (q(4)^2 + q(1:3)'*q(1:3))*eye(3) + 2*q(1:3)*q(1:3)' - 2*q(4)*skew_matrix(q(1:3));

%C_q = [ (q(4)^2 + q(1)^2 - q(2)^2 - q(3)^2)  2(q(1)q(2) + q(4)q(3))

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