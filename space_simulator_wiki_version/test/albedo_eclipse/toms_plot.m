%% TOMS Reflectivity Data Plot for 2005
close all;              % Close all figures
clear;                  % Clear workspace variables
clc;                    % Clear command window

load ga050101-051231.mat   % Load TOMS reflectivity data for the year 2005
surface(data,'EdgeColor','none')   % Plot the surface data with no edges
colormap gray           % Set colormap to grayscale
colorbar                % Display colorbar
xlim([1 288]);          % Set limit for x-axis
ylim([1 180]);          % Set limit for y-axis
fontsize = 12;          % Define font size for labels
xlabel('X-axis cells [-]', 'FontSize', fontsize)   % Add label for x-axis
ylabel('Y-axis cells [-]', 'FontSize', fontsize)   % Add label for y-axis
