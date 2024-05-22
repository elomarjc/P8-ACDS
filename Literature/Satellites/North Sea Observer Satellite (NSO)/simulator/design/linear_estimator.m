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

est_poles = [-11 -12 -13 -14 -15 -16]*0.05;

L = place(sys_w_bias.A',sys_w_bias.C',est_poles)';

lin_fil = estim(sys_w_bias,L,1:3,1:3);

% saves the states space matrices as specified in the estimator block.
A_est = lin_fil.A;
B_est = lin_fil.B;
C_est = [zeros(3) eye(3)]; % we only want the bias estimation
D_est = zeros(3,6);          % fits the direct feedthrough to the
                             % output matix

clear A B C h Sh I sys_w_bias lin_fil L est_poles
