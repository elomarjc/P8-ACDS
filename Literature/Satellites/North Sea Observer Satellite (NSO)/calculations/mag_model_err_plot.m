%**************************************************************************
% This file is used to plot the difference in a model with or
% without the motor inductance used in the magnetorquer 
% modelling.
%
% Author: Group 06gr1032
%**************************************************************************
clc
clear all
close all
format compact

%Make eps (1=true)
printEPS = 1;

sim('mag_model_err.mdl')
figure(1)
subplot(2,1,1)
plot(mag_induct_err.time,mag_induct_err.signals(1,1).values)
legend('With Coil','Without Coil','Location','NorthWest',0);
title('Magnetorquer Model Step Responces')
ylabel('Amplitude [Am^2]')
axis([0.99 1.15 202e-3 203.05e-3])
subplot(2,1,2)
plot(mag_induct_err.time,mag_induct_err.signals(1,2).values)
title('Magnetorquer Model Error')
ylabel('Amplitude difference [Am^2]')
xlabel('Simulation time [s]')
axis([0.99 1.15 0 3e-4])

if(printEPS)
    print -depsc mag_inductance_error.eps
end