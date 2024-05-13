clc
close all

% data = csvread('data/mpu9250Data.csv',7,0);

% N = size(data(:,1),1);
% x = 1:1:N(1);

% Raw
A = figure;
plot(sat_acc);
legend('Acc.x','Acc.y','Acc.z');
title('MPU9250 Accelerometer data');
ylabel('Accelerometer raw values')
xlabel('Time [s]');
saveas(A,'plots/MPU9250_acc_raw','jpg');
matlab2tikz('tikz/acc_raw.tikz', 'height', '\figureheight', 'width', '\figurewidth');

B = figure;
plot(sat_gyro);
legend('Gyro.x','Gyro.y','Gyro.z');
title('MPU9250 Gyro data');
ylabel('Gyro raw values')
% ylim([-65536/2 65536/2])
xlabel('Time [s]');
saveas(B,'plots/MPU9250_Gyro_raw','jpg');
matlab2tikz('tikz/gyro_raw.tikz', 'height', '\figureheight', 'width', '\figurewidth');

C = figure;
plot(sat_mag);
legend('Mag.x','Mag.y','Mag.z');
title('MPU9250 Magnetometer data');
ylabel('Mag raw values')
xlabel('Time [s]');
saveas(C,'plots/MPU9250_Mag_raw','jpg');
matlab2tikz('tikz/mag_raw.tikz', 'height', '\figureheight', 'width', '\figurewidth');

D = figure;
plot(sat_rpm);
legend('RW1','RW2','RW3','RW4');
title('MAXON Motor Speed');
ylabel('RPM Raw')
xlabel('Time [s]');
ylim([0 65536])
saveas(D,'plots/maxon_raw','jpg');
matlab2tikz('tikz/rw_raw.tikz', 'height', '\figureheight', 'width', '\figurewidth');


%%
E = figure;
% Acc.x PSD
lacc = length(sat_acc.data(:,1));
M1x = mean(sat_acc.data(1:round(lacc/3),1));
M2x = mean(sat_acc.data(end-round(lacc/3):end,1));

trendx = (M2x-M1x)/lacc * (1:lacc);
data_tx = sat_acc.data(:,2)-trendx';

subplot(2,2,1);
[pxx, F] = periodogram(data_tx - mean(data_tx),[],'onesided',512,100);
plot(log10(F),20*log10(pxx))
% ylim([-80 30])
title('Acc.x PSD')

%%
% Acc.y PSD
M1y = mean(sat_acc.data(1:round(lacc/3),2));
M2y = mean(sat_acc.data(1:round(lacc/3),2));

trendy = (M2y-M1y)/lacc * (1:lacc);
data_ty = sat_acc.data(:,2)-trendy';

subplot(2,2,2);
[pyy, F] = periodogram(data_ty - mean(data_ty),[],'onesided',512,100);
plot(log10(F),20*log10(pyy))
% ylim([-80 30])
title('Acc.y PSD')

% Acc.z PSD
M1z = mean(sat_acc.data(1:round(lacc/3),3));
M2z = mean(sat_acc.data(1:round(lacc/3),3));

trendz = (M2z-M1z)/lacc * (1:lacc);
data_tz = sat_acc.data(:,3)-trendz';
subplot(2,2,3);
[pzz, F] = periodogram(data_tz - mean(data_tz),[],'onesided',512,100);
plot(log10(F),20*log10(pzz))
title('Acc.z PSD')
%ylim([-80 30])

saveas(E,'plots/MPU9250_ACC_PSD','jpg');

%%
F1 = figure;
% Gyro.x PSD
lgyro = length(sat_gyro.data(:,1));
M1x = mean(sat_gyro.data(1:round(lgyro/3),1));
M2x = mean(sat_gyro.data(end-round(lgyro/3):end,1));

trendx = (M2x-M1x)/lgyro * (1:lgyro);
data_tx = sat_gyro.data(:,1)-trendx';

subplot(2,2,1);
[pxx F] = periodogram(data_tx - mean(data_tx),[],'onesided',512,100);
plot(log10(F),20*log10(pxx))
% ylim([-100 -10])
title('Gyro.x PSD')

% Gyro.y PSD
M1y = mean(sat_gyro.data(1:round(lgyro/3),2));
M2y = mean(sat_gyro.data(end-round(lgyro/3):end,2));

trendy = (M2y-M1y)/lgyro * (1:lgyro);
data_ty = sat_gyro.data(:,2)-trendy';

subplot(2,2,2);
[pyy, F] = periodogram(data_ty - mean(data_ty),[],'onesided',512,100);
plot(log10(F),20*log10(pyy))
% ylim([-100 -10])
title('Gyro.y PSD')

% Gyro.z PSD
M1z = mean(sat_gyro.data(1:round(lgyro/3),3));
M2z = mean(sat_gyro.data(end-round(lgyro/3):end,3));

trendz = (M2z-M1z)/lgyro * (1:lgyro);
data_tz = sat_gyro.data(:,3)-trendz';

subplot(2,2,3);
[pzz, F] = periodogram(data_tz - mean(data_tz),[],'onesided',512,100);
plot(log10(F),20*log10(pzz))
title('Gyro.z PSD')
%ylim([-100 -10])

saveas(F1,'plots/MPU9250_GYRO_PSD','jpg');

%%
G = figure;
% Gyro.x PSD
lmag = length(sat_mag.data(:,1));
M1x = mean(sat_mag.data(1:round(lmag/3),1));
M2x = mean(sat_mag.data(end-round(lmag/3):end,1));

trendx = (M2x-M1x)/lmag * (1:lmag);
data_tx = sat_mag.data(:,1)-trendx';

subplot(2,2,1);
[pxx, F] = periodogram(data_tx - mean(data_tx),[],'onesided',512,100);
plot(log10(F),20*log10(pxx))
% ylim([-100 -10])
title('mag.x PSD')

% Gyro.y PSD
M1y = mean(sat_mag.data(1:round(lmag/3),2));
M2y = mean(sat_mag.data(end-round(lmag/3):end,2));

trendy = (M2y-M1y)/lmag * (1:lmag);
data_ty = sat_mag.data(:,2)-trendy';

subplot(2,2,2);
[pyy, F] = periodogram(data_ty - mean(data_ty),[],'onesided',512,100);
plot(log10(F),20*log10(pyy))
% ylim([-100 -10])
title('Mag.y PSD')

% Gyro.z PSD
M1z = mean(sat_mag.data(1:round(lmag/3),3));
M2z = mean(sat_mag.data(end-round(lmag/3):end,3));

trendz = (M2z-M1z)/lmag * (1:lmag);
data_tz = sat_mag.data(:,3)-trendz';

subplot(2,2,3);
[pzz, F] = periodogram(data_tz - mean(data_tz),[],'onesided',512,100);
plot(log10(F),20*log10(pzz))
title('Mag.z PSD')
%ylim([-100 -10])

saveas(G,'plots/MPU9250_MAG_PSD','jpg');

%%
% Scale factor matrices for the sensor models
ACC_RANGE = 2;
GYRO_RANGE = 250;
MAG_RANGE = 4800;
RPM_RANGE = 16000;
DEG2RAD = pi/180;
RPM2RAD_S = 2*pi/60;

K_acc = ACC_RANGE/(2^15)                % +/- 4g int16 resolution
K_gyro = GYRO_RANGE*DEG2RAD/(2^15)              % +/- 17.4533 rad/s (1000 deg/s) int16 resolution
K_mag = MAG_RANGE/(2^15)                % +/- 4800 uT int16 resolution
K_rpm = RPM_RANGE*RPM2RAD_S/(2^16)               % + RPM uint16 resolution

%
mu_acc.x = mean(sat_acc.data(:,1)*K_acc)
var_acc.x = var(sat_acc.data(:,1)*K_acc);

mu_acc.y = mean(sat_acc.data(:,2)*K_acc)
var_acc.y = var(sat_acc.data(:,2)*K_acc);

mu_acc.z = mean(sat_acc.data(:,3)*K_acc)
var_acc.z = var(sat_acc.data(:,3)*K_acc);

mu_gyro.x = mean(sat_gyro.data(:,1)*K_gyro)
var_gyro.x = var(sat_gyro.data(:,1)*K_gyro);

mu_gyro.y = mean(sat_gyro.data(:,2)*K_gyro)
var_gyro.y = var(sat_gyro.data(:,2)*K_gyro);

mu_gyro.z = mean(sat_gyro.data(:,3)*K_gyro)
var_gyro.z = var(sat_gyro.data(:,3)*K_gyro);

mu_mag.x = mean(sat_mag.data(:,1)*K_mag)
var_mag.x = var(sat_mag.data(:,1)*K_mag);

mu_mag.y = mean(sat_mag.data(:,2)*K_mag)
var_mag.y = var(sat_mag.data(:,2)*K_mag);

mu_mag.z = mean(sat_mag.data(:,3)*K_mag)
var_mag.z = var(sat_mag.data(:,3)*K_mag);

mu_rw.a = mean(sat_rpm.data(:,1)*K_rpm)
var_rw.a = var(sat_rpm.data(:,1)*K_rpm);

mu_rw.b = mean(sat_rpm.data(:,2)*K_rpm)
var_rw.b = var(sat_rpm.data(:,2)*K_rpm);

mu_rw.c = mean(sat_rpm.data(:,3)*K_rpm)
var_rw.c = var(sat_rpm.data(:,3)*K_rpm);

mu_rw.d = mean(sat_rpm.data(:,4)*K_rpm)
var_rw.d = var(sat_rpm.data(:,4)*K_rpm);

p = [ones(1,3)*max([var_acc.x var_acc.y var_acc.z])...
     ones(1,3)*max([var_gyro.x var_gyro.y var_gyro.z])...
     ones(1,3)*max([var_mag.x var_mag.y var_mag.z])...
     ones(1,4)*max([var_rw.a var_rw.b var_rw.c var_rw.d])]'

% Sensor_gain = diag([K_acc K_gyro K_mag K_rpm]);
% P = P*Sensor_gain



