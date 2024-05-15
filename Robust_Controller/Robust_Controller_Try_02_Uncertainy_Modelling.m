close all
clc
addpath("functions\")

%% Uncertainty: There is a 5% parametric uncertainty in the inertia matrix
% Parametric uncertainty in the inertia matrix
k1 = ureal("k1",1,'Percentage',10);
k2 = ureal("k2",1,'Percentage',10);
k3 = ureal("k3",1,'Percentage',10);

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
N = Nu+K*Nx;

sys_tracking = ss(sys_ol.A-sys_ol.B*K,sys_ol.B*N,sys_ol.C,sys_ol.D)

[stabmarg,wcu] = robstab(sys_tracking);
stabmarg
stepplot(sys_tracking);

[stabmarg,wcu] = robstab(sys_tracking);
stabmarg

%% What about the integral controller?
Nr = inv([sys_ol.A,sys_ol.B;[zeros(3),eye(3)],zeros(3)])*[zeros(6,3);eye(3)];
Nx = Nr(1:6,:);
Nu = Nr(7:9,:);
N = Nu+Kp*Nx;

[stabmarg,wcu] = robstab(sys_tracking);
stabmarg

sys_int_tt = ss([sys_ol.A-sys_ol.B*Kp, -sys_ol.B*Ki;[zeros(3),eye(3),zeros(3)]],[sys_ol.B*N;-eye(3)],eye(9),zeros(9,3));
stepplot(sys_int_tt)
