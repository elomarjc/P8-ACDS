
addpath("functions\");



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

step(sys_fl) % It appears to work!

%% Integral controller
% Lets seee booiii
A_int = [A,zeros(6,3);
        zeros(3),eye(3),zeros(3)];
B_int = [B;zeros(3,3)];


Q = diag([0.05,0.05,0.05,0.001,0.001,0.001,1e-3,1e-3,1e-3]);
R = diag([1e11,1e11,1e11]);
Ko = lqr(ss(A_int,B_int,eye(9),zeros(9,3)),Q,R);

Ki= Ko(:,(end-2):end);
Kp =Ko(:,1:(end-3));
N = Nu+Kp*Nx


sys_int_tt = ss([A-B*Kp, -B*Ki;[zeros(3),eye(3),zeros(3)]],[B*N;-eye(3)],eye(9),zeros(9,3));


step(Ko*sys_int_tt*.5)



%% Testing nominal performance:

s = tf('s');

N = (eye(3)-[zeros(3),eye(3),zeros(3)]*sys_int_tt)

stepplot(N)
%%
sigma(N)
hinfnorm(N)





