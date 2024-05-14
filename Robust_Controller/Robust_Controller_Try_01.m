clear all
close all
clc

% Variables
Ip = diag([0.0088,0.0088,0.0044]);

w0= [0,0,0]';
q0 = [0,0,0]';

q_dot_0 = 1/2*skew3(w0)*q0;
w_dot_0 = inv(Ip)*(-skew3(w0)*Ip*w0);

A = [inv(Ip)*(skew3(Ip*w0)-skew3(w0)*Ip),    zeros(3,3);
                          1/2*(skew3(q0)+eye(3)),1/2*skew3(w0)];
B = [ -inv(Ip);
     zeros(3,3)];
C = eye(6,6);
D = zeros(6,3);

sys = ss(A,B,C,D)
%% Is the system controllable?
[U,d,V] = svd(ctrb(sys));
diag(d)

%Yes!

%% Try to control it

Q = diag([0.001,0.001,0.001,0.05,0.05,0.05]);
R = diag([1e8,1e8,1e8]);

K = lqr(sys,Q,R);



Nr = inv([A,B;[zeros(3),eye(3)],zeros(3)])*[zeros(6,3);eye(3)];
Nx = Nr(1:6,:);
Nu = Nr(7:9,:);
N = Nu+K*Nx;


sys_fl = ss(A-B*K,B*N,C,D);

step(K*sys_fl) % It appears to work!




function [skew]=skew3(u)
    skew = cross(repmat(u,1,3),eye(3));
end


