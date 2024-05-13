clc
close all
clear all

% Load data - that doesn't EXIST hihi ;-;
load sensor-data.mat

% figure;
% ekf_q_s.plot;% hold on;
% title('');
% xlabel('Time [s]');
% ylabel('Quaternion [-]');
% legend('');
% 
% figure
% sm_q_error.plot;% hold off;
% title('');
% ylabel('Quaternion [-]');
% xlabel('Time [s]');

%% (True) Roll Pitch and Yaw from sensor measurements
%figure(1);

% calibrate magnetometer
mag = sat_mag.Data(:,:)';
mag_calib_matrix = [0.6339   -0.0215    0.0031;
                    0         0.5742    0.0418;
                    0         0         0.5785]*1.0e-03;
mag_calib_bias = [-8.8193 130.0580 202.0582]';

for i = 1:size(mag,2)
    mag_calibrated(:,i) = mag_calib_matrix*(mag(:,i) - mag_calib_bias);
end

% calculate roll, pitch and yaw from sensor data
for i = 1:length(sat_acc.Data(:,1))
    roll_t(i) = rad2deg(atan2(sat_acc.Data(i,2), sat_acc.Data(i,3)));
end
roll = smooth(roll_t(:), 125);

for i = 1:length(sat_acc.Data(:,2))
    pitch_t(i) = rad2deg(atan2(sat_acc.Data(i,1), sqrt(sat_acc.Data(i,2)^2 + sat_acc.Data(i,3)^2)));
end
 pitch = smooth(pitch_t(:),100);

for i = 1:length(sat_mag.Data(:,1))
    yaw(i,1) = rad2deg(atan2(mag_calibrated(1,i), mag_calibrated(2,i)));
end

true = [roll pitch yaw];

ts = size(true(1:870,1),1)/10;
x1 = 0:0.1:(ts-0.1);

% figure
% plot(x1, true(1:870,:));% hold on;
% legend('Roll','Pitch','Yaw');
% xlabel('Time [s]');
% ylabel('Angle [deg]');
% title('True euler angles');

% ekf_q_s.plot
% hold off;
%% 

% Extract quaternion data from the loaded data
for i = 1:size(ekf_q_s.Data,3)
    ekf_q(1,i) = ekf_q_s.Data(4,:,i);
    ekf_q(2,i) = ekf_q_s.Data(1,:,i);
    ekf_q(3,i) = ekf_q_s.Data(2,:,i);
    ekf_q(4,i) = ekf_q_s.Data(3,:,i);
end

for i = 1:size(ekf_q,2)
    a = quat2eul(ekf_q(:,i)');
    deg_g_s_t(:,i) = rad2deg(a);
    deg_g_s(:,i) = deg_g_s_t(:,i)';
end

% Reorder Euler angles  - fusk
tmp = deg_g_s(1,:);
deg_g_s(1,:) = deg_g_s(3,:);
deg_g_s(3,:) = tmp;

% Reorder Euler angles - mere fusk
tmp = deg_g_s(1,:);
deg_g_s(1,:) = deg_g_s(2,:);
deg_g_s(2,:) = tmp;

% figure
% plot(x1, deg_g_s(1,1:870), x1, deg_g_s(2,1:870), x1, deg_g_s(3,1:870))
% legend('x','y','z');

%% Calculate estimation error
for i = 1:size(deg_g_s,2)
    est_error(:,i) = true(i,:)' - deg_g_s(:,i);
    est_error_norm(i) = sqrt(est_error(:,i)'*est_error(:,i));
end

% figure;
% plot(x1, true(1:870,1),x1, true(1:870,2),x1, true(1:870,3))
% title('true');
% 
% figure;
% plot(x1, deg_g_s(1,1:870),x1, deg_g_s(2,1:870),x1, deg_g_s(3,1:870))
% title('est');
% 
% figure;
% plot(x1, est_error(1,1:870),x1, est_error(2,1:870),x1, est_error(3,1:870))
% title('est error');

% figure;
% plot(est_error_norm)
% title('error norm');

ts = size(est_error_norm(1:855),2)/10;
x = 0:0.1:(ts-0.1);

figure
plot(x,est_error_norm(1:855))
%legend('Estimation error');
title('');
xlabel('Time [s]');
ylabel('Angular Error [deg]');
grid on;
grid minor;

matlab2tikz('tikz/attitude_error_est_norm.tikz', 'height', '\figureheight', 'width', '\figurewidth');

%figure;
%plot(x,true_euler_data(1:855))

%figure;
%plot(x,deg_g_s(1:855))



