clc

s = tf('s');

err = sind(1/2); % Half a degree err
Mpmax = 0.25; % 25% overshoot
ts_max = 145;

%% Make a plot with the required constraints
% First order ss error
w_s = s/err;

% Overshoot and settling time
zeta = -log(Mpmax)/sqrt(pi^2+(log10(Mpmax))^2)
wn = 3/ts_max/zeta
w_Mp = wn*sqrt(1-zeta^2)

W_log10_lim = [-2.5;-1]; 
W = logspace(W_log10_lim(1),W_log10_lim(2),100);



f1 = figure(1);
clf(f1);
hold on
grid on
title("Performance weight requirements");
bodemag(w_s,W,'--');
plot([10.^(W_log10_lim),[w_Mp;w_Mp]],[repmat(log10(1+Mpmax)*20,2,1),[-60;60]],'k--','DisplayName',"T4.b");
plot([10.^(W_log10_lim),[wn;wn]],[[0;0],[-60;60]],'--','Color',[1,0.5,0],'DisplayName',"T4.a");
xlim([10.^(W_log10_lim)])
ylim([-20,20])
legend

%w_1 = 0.001;
w_1 = 0.005;
dist = 0.0025;
w_s = (s/err)/(s/(w_1+dist)+1)^2*(s/w_1+1)^1;


bodemag(w_s,'r')
%% Checking the requirement:
f2 = figure(2);
clf(f2);
bodemag([1,0,0]*N_ki*[1;0;0]);
hold on
bodemag(w_s)

%%bodemag(s/err/(s/0.00863+1)*(s/0.044+1)/(s/0.1+1))