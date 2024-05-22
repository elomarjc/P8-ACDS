% Design of an attitude controller using LQR control theory.
% clear all

ang_robust_control;
operating_point_constants;

analyse = 1;

A_sat=[inv(I)*Sh];
B_sat=-inv(I);
C_sat=eye(3);
D_sat=zeros(3);

K = D_ang_c;

SYS=ss(A_sat,B_sat,C_sat,D_sat);
CL=feedback(K*SYS,eye(3),+1);

A_att = [zeros(3) 0.5*eye(3); zeros(3) CL.a];
B_att = [zeros(3); CL.b];
C_att = [eye(3) zeros(3)];%;zeros(3) CL.c];
D_att = zeros(3);

%introducing integral state and augmenting the system
A_int = [A_att zeros(6,3); C_att zeros(3)];
B_int = [B_att; zeros(3)];
C_int = eye(9);%[C_att zeros(3)];
D_int = zeros(9,3);


%weight matrix punishing the state in the system
wmax=0.0008; 
qmax=0.00552; 
integral = 0.000045; 
Q=blkdiag(qmax*eye(3),wmax*eye(3),integral*eye(3));

%weight matrix punishing the input to the system
R=1*eye(3);

%optimal lqr controller
L=-lqr(A_int,B_int,Q,R);

if analyse
L
SYS_int=ss(A_int,B_int,C_int,D_int);
CL_int=feedback(L*SYS_int,eye(3),+1);
damp(CL_int)
end
% State space model of controller
A_att_c = zeros(9,9);
B_att_c = zeros(9,9);
C_att_c = zeros(3,9); 
D_att_c = L;

if(analyse==0)
clear A_sat B_sat C_sat D_sat K SYS CL A_att B_att C_att D_att wmax qmax Q R L SYS_att CL_att analyse I Sh h;
end