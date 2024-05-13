% Clear command window
clc

% Close all figures
close all

% Clear workspace variables
clear all

% Load data from the file 'ekf_sim_data-160620.mat'
load ekf_sim_data-160620.mat

% Plot quaternion components from the EKF simulation
subplot(2,1,1)
ekf_q_s.plot
legend('q1','q2','q3','q4');
title('');
ylabel('[-]');
xlabel('');

% Plot estimation error from the EKF simulation
subplot(2,1,2)
ekf_est_error.plot
legend('Estimation error');
title('');
ylabel('Error [deg]');
xlabel('Time [s]');

% Convert the plot to TikZ format for use in LaTeX
matlab2tikz('test/satball/tikz/ekf_sim_q_og_est_error_slides.tikz', 'height', '\figureheight', 'width', '\figurewidth');
