%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build a state space representation of a kalman
% estimator which estimates the angular velocity
% caused by the slow changing disturbance torques.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

operating_point_constants;

% System matrix, input matrix and output matrix including angular
% velocity bias
A = [inv(I)*Sh inv(I);zeros(3,6)];
B = [-inv(I);zeros(3)];
C = [eye(3) zeros(3)];

% state space model of the system with bias
sys_w_bias = ss(A,B,C,0);

% state space model expanded with the stocastic process noise
sys_w_expanded = ss(sys_w_bias.A,[sys_w_bias.B eye(6)],sys_w_bias.C,0);

% process noise covariance
cov_omega = 5e-15; % 1e-11
cov_bias  = 5e-15; % 1e-12

% measurement noise covariance
cov_meas_omega = 1e-12; % 3e-12

% process and measurement noise covariane matrix
cov_proc = [cov_omega*eye(3),zeros(3);zeros(3),cov_bias*eye(3)];
cov_meas = [cov_meas_omega*eye(3)];

% calculation of the kalman filter
[kal_fil, L, P] = kalman(sys_w_expanded,cov_proc,cov_meas);

% saves the states space matrices as specified in the estimator block.
A_kalman = kal_fil.A;
B_kalman = kal_fil.B;
C_kalman = [zeros(3) eye(3)]; % we only want the bias estimation
D_kalman = zeros(3,6);          % fits the direct feedthrough to the
                                % output matix
clear A B C h Sh I sys_w_bias sys_w_expanded cov_meas_omega ...
    cov_proc cov_meas kal_fil L P cov_omega cov_bias