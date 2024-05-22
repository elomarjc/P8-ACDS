% Design of an attitude controller using LQR control theory.

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
C_att = [eye(3) zeros(3);zeros(3) CL.c];
D_att = zeros(6,3);

%weight matrix punishing the state in the system
wmax=0.01;   
qmax=0.0085;
Q=[qmax*eye(3) zeros(3); zeros(3) wmax*eye(3)];

%weight matrix punishing the input to the system
R=0.5*eye(3);

%optimal lqr controller
L=-lqr(A_att,B_att,Q,R);

if(analyse==1)
SYS_att=ss(A_att,B_att,C_att,D_att);
CL_att=feedback(L*SYS_att,eye(3),+1);
damp(CL_att)
end
% State space model of controller
A_att_c = zeros(6,6);
B_att_c = zeros(6,6);
C_att_c = zeros(3,6); 
D_att_c = L;

clear A_sat B_sat C_sat D_sat K SYS CL A_att B_att C_att D_att wmax qmax Q R L SYS_att CL_att analyse I Sh h;
