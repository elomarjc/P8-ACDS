function [u,q_error,s] = sm_gomx3(state,ref,q_s_c,q_th_c,K,delta_0,e,max_disturbance,Jsat,Jrw)
%% Quaternion Sliding Mode control for GOMX3
% state - state vector [q_i_s; omega_i_s; omega_m] (omega_m is the motor
% velocity [omega1 omega2 omega3 omega4]'
% ref - reference vector [q_ref; omega_i_ref; omega_i_ref_dot]
% q_s_c - is the rotation from satellite body to control frame
% q_th_c - is the rotation from satellite body to the tetrahedron config frame
% K - is the gain on the switching surface
% delta_0 - is the switching gain
% e - is the saturation gain
% max_disturbance - Upperbound on the disturbance
% Jsat - Satellite inertia [Ixx Iyy Izz]'
% Jrw - Reaction wheel inertia Izz

%% Define system variables
Js = diag(Jsat);
Td = max_disturbance;                               % Upperbound on the disturbance torques
Js_inv_upper = 1/min([Js(1,1),Js(2,2),Js(3,3)]);    % Upperbound of the inverse moment of inertia
v = zeros(3,1);
alpha = deg2rad(60);
beta = deg2rad(19.47);

%% Tetrahedron distribution matrix [m1 m2 m3 m4]
P_w_th = [cos(beta)    -cos(beta)*cos(alpha)  -cos(beta)*cos(alpha)    0;
          0            cos(beta)*cos(alpha/2) -cos(beta)*cos(alpha/2)  0;
          -sin(beta)   -sin(beta)             -sin(beta)               1];
      
%% Define the states
q_i_s = state(1:4);                                 % quaternion attitude of the satellite body in the inertial reference frame
omega_i_s = state(5:7);                             % Angular velocity of the satellite body relative to the inertial reference frame
omega_m = state(8:11);
%% Define the reference
q_i_ref = ref(1:4);                                 % quaternion target reference frame given in the inertial reference frame (fx. orbit reference fram)
omega_i_ref = ref(5:7);                             % The angular velocity of the target reference frame
omega_i_ref_dot = ref(8:10);                        % The angular acceleration of the target reference frame

%% Generate the quaternion error between the satellite body and the target
q_i_s = qunit(q_i_s);                               % Check form unity
q_s_e = qmult(q_i_ref,qinv(q_i_s));                 % Error quaternion

omega_s_e = omega_i_ref - omega_i_s;                % Find the angular velocity error

% %% If the coriolis moments should be included
% omega_i_tetra = qRot(omega_i_s,q_s_tetra);
% q_th_c = qmult(qinv(q_s_tetra),q_s_c);
% 
% omega_w1 = qRot(omega_i_tetra,qinv(A1_q)) + [0 0 omega_m(1)]';
% omega_w2 = qRot(omega_i_tetra,qinv(A2_q)) + [0 0 omega_m(2)]';
% omega_w3 = qRot(omega_i_tetra,qinv(A3_q)) + [0 0 omega_m(3)]';
% omega_w4 = qRot(omega_i_tetra,qinv(A4_q)) + [0 0 omega_m(4)]';
% 
% if coriolis_on
%     tau_th_p = (qRot(-cross(omega_w1,Jrw*omega_w1),A1_q)  +...
%         qRot(-cross(omega_w2,Jrw*omega_w2),A2_q)) +...
%         qRot(-cross(omega_w3,Jrw*omega_w3),A3_q)  +...
%         qRot(-cross(omega_w4,Jrw*omega_w4),A4_q);
% else
%     tau_th_p = (qRot(-cross(omega_w1,Jrw*[0 0 omega_m(1)]'),A1_q)  +...
%         qRot(-cross(omega_w2,Jrw*[0 0 omega_m(2)]'),A2_q)) +...
%         qRot(-cross(omega_w3,Jrw*[0 0 omega_m(3)]'),A3_q)  +...
%         qRot(-cross(omega_w4,Jrw*[0 0 omega_m(4)]'),A4_q);
% end
%% Find the quaternion velocity
qe_s_dot = [qTrans(q_s_e)*omega_s_e;
            -0.5 * q_s_e(1:3)'*omega_s_e];

%% Define the Sliding surface s
s = qe_s_dot(1:3) + K*q_s_e(1:3);

%% Define the size of the switching gain g0
g0 = 0.5*Td*Js_inv_upper + delta_0;

%% Saturation function
% e is the slope of the saturation function
% u0 is the switching gain
if norm(s) > e
    v = g0*s/norm(s);
elseif norm(s) <= e
    v = g0*s/e;
end

%% Define the quaternion error dynamics fb independent of the control u
omega_i_c = qRot(omega_i_s,q_s_c);

%% Coriolis moment from the reaction wheels
h_wheels = Jrw*omega_m;
tau_c_p = cross(omega_i_c,qRot(P_w_th*h_wheels,q_th_c));
 
fb = qTrans(q_s_e)*( omega_i_ref_dot -...
    qRot(Js^-1 * (-cross(omega_i_c,Js*omega_i_c) - tau_c_p),qinv(q_s_c)))...
    - 0.25*q_s_e(1:3)*(omega_s_e'*omega_s_e);

%% Define the control law u which drives the system to the sliding surface
u = -(-qRot(Js*qRot(qTrans(q_s_e)'*(fb + K*qe_s_dot(1:3) +...
    v),q_s_c),qinv(q_s_c)));

% Output the quaternion error for reference
q_error = q_s_e;
end

% function [ p ] = qmult( q1,q2 )
% %UNTITLED2 Summary of this function goes here
% %   Detailed explanation goes here
% p = [q1(4)*q2(1:3) + q2(4)*q1(1:3) + cross(q1(1:3),q2(1:3));
%     q1(4)*q2(4) - dot(q1(1:3),q2(1:3))];
% 
% end

% function [ q_inv ] = qinv( q )
% %UNTITLED Summary of this function goes here
% %   Detailed explanation goes here
% q_inv = [-q(1) -q(2) -q(3) q(4)]';
% end
% 
% function [ v ] = qRot( u,q )
% %UNTITLED4 Summary of this function goes here
% %   Detailed explanation goes here
%     v_hat = qmult(qmult(qinv(q),[u;0]),q);
%     v = v_hat(1:3);
% 
% end
% 
% function [ T ] = qTrans( q )
% %UNTITLED5 Summary of this function goes here
% %   Detailed explanation goes here
% T = 0.5*[q(4) -q(3)  q(2);
%         q(3)  q(4) -q(1);
%         -q(2)  q(1)  q(4)];
% end

% function [ v ] = qunit( q )
% %UNTITLED3 Summary of this function goes here
% %   Detailed explanation goes here
% 
% v = q./norm(q);
% end

