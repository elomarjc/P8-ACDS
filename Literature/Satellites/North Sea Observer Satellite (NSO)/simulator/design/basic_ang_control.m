% Design of state feedback control for the angular velocity of the
% satellite, the operating point of the angular velocity is zero and the
% for the angular momentum the bias point of 4.5e-5 is chosen as operating
% point. Since the operating point of the angular velocity is zero the
% change in angular momentum is left out of the system equation, thus A and
% B is given as

% operating point
operating_point_constants;

% System description
A_basic_ang=[inv(I)*Sh];
B_basic_ang=-inv(I);

size_A_basic_ang = size(A_basic_ang);

CO_basic_ang=ctrb(A_basic_ang,B_basic_ang);
if rank(CO_basic_ang) < size_A_basic_ang(1)
  error(strcat('The system is not controlable. States: ',size_A_basic_ang(1),' rank of controllability matrix: ',rank(CO)))
end

% Desired poles
P_basic_ang=[-1.254 -1.255 -1.256]'; % placing the poles

K_basic_ang=-place(A_basic_ang,B_basic_ang,P_basic_ang);


% State space model of controller
A_ang_c = zeros(3,3);
B_ang_c = zeros(3,3);
C_ang_c = zeros(3,3); 
D_ang_c = K_basic_ang;

clear K_basic_ang P_basic_ang CO_basic_ang size_A_basic_ang A_basic_ang B_basic_ang I Sh h