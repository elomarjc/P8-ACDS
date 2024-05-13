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

% Estimeret roll, pitch og yaw
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

% fusk
tmp = deg_g_s(1,:);
deg_g_s(1,:) = deg_g_s(3,:);
deg_g_s(3,:) = tmp;

% mere fusk
tmp = deg_g_s(1,:);
deg_g_s(1,:) = deg_g_s(2,:);
deg_g_s(2,:) = tmp;

ts = size(deg_g_s(1:730),2)/10;
x = 0:0.1:(ts-0.1);

figure;
plot(x',deg_g_s(:,1:730)'); hold on;
plot([0 73], [-80 -80], '--k'); hold off;
legend('Roll','Pitch','Yaw', 'Yaw Reference');
ylabel('Angle [deg]');
xlabel('Time [s]');
grid on;
grid minor;

matlab2tikz('tikz/80_deg_yaw_step_sm.tikz', 'height', '\figureheight', 'width', '\figurewidth');





