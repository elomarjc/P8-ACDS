% operating point
h_bias = 4.139e-6*1078.6/2;
h=[h_bias h_bias h_bias];
Sh=[0 -h(3) h(2);h(3) 0 -h(1);-h(2) h(1) 0];
I=diag([42.3 42.3 28.4]/1e3); % inertia matrix for deployed situation

% System description
A=[inv(I)*Sh];
B=-inv(I);
C=eye(3);
sys = ss(A,B,C,0);

M = 1;
A = 0.01;
wb = 0.01;

w_ks = 10;
M_ks = 1;
A_ks = 0;

w_wt = 0.5;
M_wt = 1;
A_wt = 1;

W1 = tf([1/M wb],[1 wb*A]);
W3 = tf([1 w_wt/M_wt],[A_wt w_wt]);
W2 = tf([1 0],[1 w_ks]);

% $$$ w0=1;     % desired closed-loop bandwidth
% $$$ A=1/100;  % desired disturbance attenuation inside bandwidth
% $$$ M=2 ;      % desired bound on hinfnorm(S) & hinfnorm(T)
% $$$ W1=tf([1/M w0],[1 w0*A]); % Sensitivity weight
% $$$ W2=0;                % Empty control weight
% $$$ W3=(s+w0/M)/(A*s+w0); % Complementary sensitivity weight


[K1,CL,Gam,Info] = mixsyn(sys,W1,[],W3);
K = K1;
% State space model of controller
A_ang_c = K.A;
B_ang_c = K.B;
C_ang_c = K.C;
D_ang_c = K.D;

L=sys*K;  % loop transfer function
S=inv(eye(3)+L); % Sensitivity
T=eye(3)-S;      % complementary sensitivity
sigma(Gam/W1,'g-.',Gam*sys/ss(W2),'r-.',S,'g',T,'r');%,L,'b',Gam/W2,'b-.');
