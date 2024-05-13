clc                 % Clear command window
clear all           % Clear all variables from workspace
close all           % Close all figures

Ts = 1;             % Sampling time
Ts_motor = 1/4;     % Motor sampling time
alpha = deg2rad(60);    % Angle in radians
beta = deg2rad(19.47);  % Angle in radians

% Transformation matrix from wheel to tetrahedron
P_w_th = [cos(beta)    -cos(beta)*cos(alpha)  -cos(beta)*cos(alpha)    0;
          0            cos(beta)*cos(alpha/2) -cos(beta)*cos(alpha/2)  0;
          -sin(beta)   -sin(beta)             -sin(beta)               1];

P_th_w = P_w_th'*(P_w_th*P_w_th')^-1;   % Inverse transformation matrix

% Moment of inertia matrix
Jb = [6.714e2 -9.903e-1 7.125e-1;
    -9.903e-1 7.044e2 9.894e-1;
    7.125e-1 9.894e-1 7.267e2] * 1.0e-6;

[R_s_c,Js] = eig(Jb);  % Eigenvalues and eigenvectors

Jsx = Js(1,1);  % Eigenvalue along x-axis
Jsy = Js(2,2);  % Eigenvalue along y-axis
Jsz = Js(3,3);  % Eigenvalue along z-axis

%% Motor Specifications
Kt = 1.81e-3;   % Motor torque constant
Ke = 1.81e-3;   % Motor back EMF constant

Rmo = 4.44;     % Motor resistance
bm = 6.3895e-08;    % Motor viscous damping coefficient

Lambda_m = (Kt*Ke/Rmo + bm);    % Motor damping coefficient
Gamma_m = (Kt/Rmo);             % Motor torque constant
Jm = 0.3811e-6;                  % Motor inertia

%% Reference frames
q_s_c = A2q(R_s_c);             % Convert rotation matrix to quaternion
q_s_tetra = [0 0 1 0]';         % Quaternion for tetrahedron frame
q_th_c = qmult(qinv(q_s_tetra),q_s_c);   % Quaternion for tetrahedron to wheel frame
R_th_c = quat2rotm([q_th_c(4) q_th_c(1) q_th_c(2) q_th_c(3)])'; % Rotation matrix

Wb_o = [0 -0.0012 0]';  % Angular velocity vector in the body frame
Wc_o = qRot(Wb_o,q_s_c);   % Angular velocity vector in the wheel frame

%% Motor model
Am = -Jm^-1 * Lambda_m; % Motor state matrix
Bm = Jm^-1 * Gamma_m;   % Motor input matrix
Cm = 1;                 % Motor output matrix
Dm = 0;                 % Motor feedthrough matrix

ss_motor = ss(Am,Bm,Cm,Dm); % State-space representation of motor model
[b,a] = ss2tf(Am,Bm,Cm,Dm); % Convert state-space model to transfer function
tf_motor = tf(b,a);         % Transfer function of motor
D_motor = c2d(tf_motor,Ts_motor,'Tustin');    % Discretize motor transfer function
c = pid(0.0018,0.004);      % Define PID controller
cl = feedback(c*tf_motor,1);    % Closed-loop transfer function
d_cl = c2d(cl,Ts_motor,'zoh');  % Discretize closed-loop transfer function

[y,t] = step(d_cl*1000,5);     % Step response of closed-loop system

% Plot motor response
figure
stairs(t,y)
hold on 
plot([0 5],[1000 1000],'--')
grid on
title('Motor PI Regulator')
xlabel('Time [s]')
ylabel('Angular velocity [rad/s]')
legend('Motor speed','Speed Reference');
ylim([0 1500]);

pole_motor_of_cl = pole(cl)   % Calculate poles of closed-loop system

%% Calculation for Angular Acceleration
Aws = -Js^-1 * (skew_matrix(Js*Wc_o) - skew_matrix(Wc_o)*Js);   % Angular acceleration matrix
Aq = 0.5*skew_matrix(Wb_o);    % Quaternion derivative matrix

A = [Aq 0.5*eye(3);
    zeros(3)  Aws];    % System matrix

B = [zeros(3);-Js^-1];  % Input matrix
C = eye(6);             % Output matrix
D = zeros(6,3);         % Feedthrough matrix

ss_alpha = ss(A,B,C,D); % State-space representation for angular acceleration
D_alpha = c2d(ss_alpha,Ts,'Tustin');    % Discretize state-space system
rank_of_ss_alpha = rank(ctrb(ss_alpha))    % Rank of controllability matrix

%% LQR Controller Design
Q = [ones(1,3)*deg2rad(0.5)^-2 ones(1,3)*deg2rad(1)^-2];   % State error weights
Q = diag(Q);
R = (eye(3)*(3e-6)^(-2));   % Control input weights
[K_lqr,S,e] = dlqr(D_alpha.a,D_alpha.b,Q,R);   % LQR controller gain

save('lqr_gain.mat','K_lqr');  % Save LQR gain

alpha_cl = ss_alpha.a-ss_alpha.b*K_lqr;    % Closed-loop system matrix

sys_cl = ss(alpha_cl,ss_alpha.b,ss_alpha.c,ss_alpha.d); % Closed-loop state-space system
pole_sys = abs(pole(sys_cl))         % Magnitude of poles of closed-loop system
pole_motor_of_abs_cl = abs(pole(cl)) % Magnitude of poles of closed-loop motor system
