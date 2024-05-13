clc
clear all
close all
q = [0 0 0 1]';
omega = [0 0 pi]';
dt = 0.001;

for i = 1:1000
    d_q_hat = qTrans(q)*omega*dt;
    d_q_4 = -1/2 * q(1:3)'*omega*dt;
    
    q = q + [d_q_hat;d_q_4];
%     q = q/norm(q);
end

