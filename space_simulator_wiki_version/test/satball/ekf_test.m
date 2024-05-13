clc
close all
clear all

% Initial quaternion for the satellite attitude estimation
q_s_c = qunit([0.0218   -0.0064   -0.0153    0.9996]); % Corrected quaternion
q_s_tetra = [0 0 -1 0]; % Quaternion for tetrahedral arrangement

% Initial states and biases
q_s_0 = [0 0 0 1]; % Initial quaternion
omega_0 = [0 0 0]; % Initial angular velocity
acc_bias_0 = [0 0 0]; % Initial accelerometer bias
gyro_bias_0 = [0 0 0]; % Initial gyroscope bias
mag_bias_0 = [0 0 0]; % Initial magnetometer bias
omega_rw = [800 800 800 800]; % Gyroscope random walk

% Initial state vector
x0 = [q_s_0 omega_0 omega_rw]';

% Control input
u = ones(4,1)*1.5; % Constant control input

% Sensor measurements
acc = [0 0 1]; % Accelerometer measurement
gyro = [0 0 0]; % Gyroscope measurement
mag = [1 0 0]; % Magnetometer measurement
rw_omega = [800 800 800 800]*0; % Gyroscope white noise

% Measurement vector
z_meas = [acc gyro mag rw_omega]';

% Time step
dt = 1;

% EKF iteration
for i = 1:1
    [x_hat,phi,H] = satball_ekf(z_meas,u,x0,dt,q_s_c',q_s_tetra');
    
    % Observability analysis
    Ob = rank(obsv(phi,H)); % Rank of observability matrix
    unob(i) = length(phi) - Ob; % Unobservable states count
end
