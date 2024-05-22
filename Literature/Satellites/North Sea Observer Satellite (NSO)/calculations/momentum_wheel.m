%**************************************************************************
% This file is used to plot some illustrations used in the momentum wheel 
% modelling. Using the motor constants the difference in a model with or
% without the motor inductance is plotted. Furthermore, a PI controller
% with maximum gain is calculated.
%
% Author: Group 06gr1032
%**************************************************************************
clc
clear all
close all
format compact

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ts = 409.67e-6; %1/(fs)=1/2.441e3 
s = tf('s');
printEPS = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Motor model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Motor constants
R=8; %8ohm
L=0.07e-3; %0.07mH
Kt=2.67E-3; %2.67mNm/A
Ks=1/373.85; %1/(3570/60*(2*pi)) %rpm/V -> V/(rad/s)
J=6.7E-9+4.1332e-6; %0.067g/cm^2 + 4.1332e-006 flywheel
no_load_current = 16e-3;
no_load_speed = (10300/60*2*pi); %rpm -> rad/s
b=(no_load_current*Kt)/no_load_speed; %calculate from a steady state.. 
%% J*domega/dt = Kt*i_mw - b*omega = 0 %used to be: 7E-9; 

sys = feedback(tf([Kt],[R*J b*R]),Ks);
sys_induct = feedback(tf([Kt],[J*L R*J+b*L b*R]),Ks);
[out t] = step(sys);
t = 0:Ts:t(end);
[out t] = step(sys,t);
[out_induct t_induct] = step(sys_induct,t);

%%%Model error without inductance%%%
figure(1)
hold on
subplot(2,1,1)
plot(t(2:700:end),out(2:700:end),'bx',t_induct,out_induct,'r')
legend('Without Coil','With Coil','Location','NorthWest',0);
title('Model Step Responces')
ylabel('Amplitude [rad/s]')
subplot(2,1,2)
plot(t,out-out_induct)
title('Model Error')
ylabel('Amplitude difference [rad/s]')
xlabel('Simulation time [s]')
hold off
if(printEPS)
    print -depsc momentum_model_diff.eps
end

%%%%%%%% Step with P controller esuring clp unit gain%%%%%%%%
sysc_unit = tf([0 1/dcgain(sys) 0],[1 0]); %P with clp unit gain
clp_unit = feedback(series(sys,sysc_unit),1);
clp_induct_unit = feedback(series(sys_induct,sysc_unit),1);

[out_clp_unit t_clp_unit] = step(clp_unit);
t_clp_unit = 0:Ts:t_clp_unit(end);
[out_clp_unit t_clp_unit] = step(clp_unit,t_clp_unit);
[out_induct_clp_unit t_induct_clp_unit] = step(clp_induct_unit,t_clp_unit);

figure(2)
hold on
subplot(2,1,1)
step(clp_unit,'bx',clp_induct_unit,'r')
subplot(2,1,2)
title('Model Error')
plot(t_clp_unit ,out_clp_unit-out_induct_clp_unit)
hold off

%%%%%%%% Bode of P controller %%%%%%%%
figure(3)
hold on
%set(gca,'Fontsize',12);
bode(clp_unit, logspace(-2,2.8845397,1000)) %max 122hz = 766rad/s
title('Momentum Wheel Unit Gain Closed Loop Response');
hold off
if(printEPS)
    print -depsc momentum_wheel_cloop.eps
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HW-Controller model %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gain_missing = 20*log10(evalfr(clp_unit,122*2*pi)); %122 is max bandwidth
disp(strcat('Missing: ',num2str(gain_missing),' [db]'))
Kp = 10^(-gain_missing/20)*(1/dcgain(sys));
Kp = 0.4; %on account of actuator limitations
disp(strcat('Controller P: ',num2str(Kp)))
Ki = Kp/10;
Kd = 0; %0
sysc=tf([Kd Kp Ki],[1 0]); %PID

% Discretisize
[dencz,numcz]=c2dm([1 0],[Kd Kp Ki],Ts,'tustin'); 

%%%%%%%%%%%% Motor plot %%%%%%%%%%%%%%%%%%%%%%%%% 
%rlocus(sys)
figure(4)
hold on
subplot(2,2,1)
bode(sys, logspace(-2,3.9,1000))
title('Motor W_m/V_m');
[num,den]=tfdata(sys);
[numz,denz] = c2dm(num,den,Ts,'zoh');

%%%%%%%%%%%% Controller plot %%%%%%%%%%%%%%%%%%%%%
subplot(2,2,2)
bode(sysc, logspace(-2,3.6,1000))
title('Hardware controller');
%%%%%%%%%%%% Forward coupling %%%%%%%%%%%%%%%%%%%%
subplot(2,2,3)
bode(series(sysc,sys), logspace(-2,2.8845397,1000))
title('Harware Controller * Motor');
%subplot(2,2,4)
hold off

%%%%%%%%%% Step with PI controller %%%%%%%%%%%%%%%
clp = feedback(series(sysc,sys),1);
clp_induct = feedback(series(sysc,sys_induct),1);

[out_clp t_clp] = step(clp);
t_clp = 0:Ts:t_clp(end);
[out_clp t_clp] = step(clp,t_clp);
[out_induct_clp t_induct_clp] = step(clp_induct,t_clp);

figure(5)
hold on
subplot(2,1,1)
step(clp,'bx',clp_induct,'r')
subplot(2,1,2)
title('Model Error')
plot(t_clp ,out_clp-out_induct_clp)
hold off

%%%%%%%% Bode of PI controller %%%%%%%%
figure(6)
hold on
%set(gca,'Fontsize',12);
bode(clp, logspace(-2,2.8845397,1000))
title('Momentum Wheel Closed Loop Response');
hold off
if(printEPS)
    print -depsc momentum_wheel_cloop_maxp.eps
end

%%%%%%%%%%% Motor torque test %%%%%%%%
m_bias = round(5000/60*(2*pi));%5000rpm
m_max = round(10000/60*(2*pi));%10000rpm
u = [m_bias:1:m_max m_max-1:-1:m_bias];
u = [u u(2:end) u(2:end)];
time = 10/length(u):10/length(u):10;
omega_mw = lsim(clp,u,time);
% ydif = diff(y)*Ts;
% ydif = [ydif' ydif(end)]';
% torque = J*ydif;
torque = J*lsim(clp*s,u,time);

plotstart =12;
figure(7)
subplot(2,1,1)
hold on
plot(time(plotstart:30:end),u(plotstart:30:end),'bx')
plot(time(plotstart:end),omega_mw(plotstart:end),'r')
title('Angular Velocity and Reference')
ylabel('\omega_{mw} [rad/sec]');
legend('\omega_{ref}','\omega_{mw}','Location','NorthWest',0);
subplot(2,1,2)
hold on
plot(time(plotstart:end),torque(plotstart:end))
title('Torque Produced')
ylabel('N_{mw} [Nm]');
xlabel('Simulation time [s]');
hold off
if(printEPS)
    print -depsc momentum_torque.eps
end