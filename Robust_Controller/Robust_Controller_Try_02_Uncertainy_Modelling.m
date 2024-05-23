close all
clc
addpath("functions\")

%% Uncertainty: There is a 5% parametric uncertainty in the inertia matrix
% Parametric uncertainty in the inertia matrix
k1 = ureal("k1",1,'Percentage',5);
k2 = ureal("k2",1,'Percentage',5);
k3 = ureal("k3",1,'Percentage',5);

% Neglected dynamics uncertainty caused by the linearization
w1 = ureal("w1",0,'Plusminus',0.011);
w2 = ureal("w2",0,'Plusminus',0.011);
w3 = ureal("w3",0,'Plusminus',0.011);


q1 = ureal("q1",0,'Plusminus',0.0436194);
q2 = ureal("q2",0,'Plusminus',0.0436194);
q3 = ureal("q3",0,'Plusminus',0.0436194);

uIp = diag([0.0088*k1,0.0088*k2,0.0044*k3]);

w0= [w1,w2,w3]';
q0 = [q1,q2,q3]';

q_dot_0 = 1/2*skew3(w0)*q0;
w_dot_0 = inv(uIp)*(-skew3(w0)*uIp*w0);

A = [inv(uIp)*(skew3(uIp*w0)-skew3(w0)*uIp),    zeros(3,3);
                          1/2*(skew3(q0)+eye(3)),1/2*skew3(w0)];
B = [ -inv(uIp);
     zeros(3,3)];
C = eye(6,6);
D = zeros(6,3);

sys_ol = ss(A,B,C,D);

%% Is the system robust stable?
% Using K from previous file

%margins = allmargin(sys_ol*K);

sys_cl = feedback(sys_ol*K,eye(6));

Nr = inv([sys_ol.A,sys_ol.B;[zeros(3),eye(3)],zeros(3)])*[zeros(6,3);eye(3)];
Nx = Nr(1:6,:);
Nu = Nr(7:9,:);
N_k = Nu+K*Nx;

sys_tracking = ss(sys_ol.A-sys_ol.B*K,sys_ol.B*N_k,sys_ol.C,sys_ol.D)

[stabmarg,wcu] = robstab(sys_tracking);
stabmarg
stepplot(sys_tracking);



%% What about the integral controller?
Nr = inv([sys_ol.A,sys_ol.B;[zeros(3),eye(3)],zeros(3)])*[zeros(6,3);eye(3)];
Nx = Nr(1:6,:);
Nu = Nr(7:9,:);
N_k = Nu+Kp*Nx;


sys_int_tt = ss([sys_ol.A-sys_ol.B*Kp, -sys_ol.B*Ki;[zeros(3),eye(3),zeros(3)]],[sys_ol.B*N_k;-eye(3)],eye(9),zeros(9,3));


%[stabmarg,wcu] = robstab(sys_int_tt);
%stabmarg

%%
stepplot(sys_int_tt)

%% Testing nominal performance:
s = tf('s')
%Performance weight: Step of 60 deg:
% +-5deg of 1st order -> err < 12%
err = 0.12;
% must contain 0.0257 rad/s 
w_in = 0.0257;
% must have an overshoot smaller than 15 deg -> Mp < 25% to be checked -> damping ratio >0.4 -> 1.36 (1/(2*psi*sqrt(1-psi^2)))
max_bound = 1.36;
% w_p as the S bound
wp = (err*s)/(err*s/max_bound+1);
bodemag(wp);
hold on
grid on
plot([w_in,w_in],[-60,60],'--k');
wp = blkdiag(wp,wp,wp);
%% Checking non-integral
sys_tracking_q = [zeros(3),eye(3)]*sys_tracking;

N_k = eye(3)-sys_tracking_q;
out = 1/wp*N_k;

%% Checking integral
sys_int_tt_quat = [zeros(3),eye(3),zeros(3)]*sys_int_tt;

s = tf('s');

N_ki = (eye(3)-sys_int_tt_quat);
%bodemag(N_k,N_ki,wp)
legend("K","Ki","wp");
figure
out = 1/wp*N_ki;
plot_performance(out)
%%
sigma(N_k)
hinfnorm(N_k)


function []= plot_performance(out)
    subplot(131);
    bodemag([1,0,0]*out*[1,0,0]');
    hold on
    plot([0.0001,10000],[0,0],'--k');
    subplot(132);
    bodemag([0,1,0]*out*[0,1,0]');
    hold on
    plot([0.0001,10000],[0,0],'--k');
    subplot(133);
    bodemag([0,0,1]*out*[0,0,1]');
    hold on
    plot([0.0001,10000],[0,0],'--k');
end
