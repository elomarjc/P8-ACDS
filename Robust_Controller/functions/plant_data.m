function [sys] = plant_data()
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

sys = ss(A,B,C,D);



end