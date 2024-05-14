close all
clc
addpath("functions\")

%% Uncertainty: There is a 5% parametric uncertainty in the inertia matrix
k = ureal("k",1,'Range',[0.95,1.05]);
uIp = diag([0.0088,0.0088,0.0044])*k;

w0= [0,0,0]';
q0 = [0,0,0]';

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

margins = allmargin(sys_ol*K);

sys_cl = feedback(sys_ol*K,eye(6));
[stabmarg,wcu] = robstab(sys_cl);
stabmarg
