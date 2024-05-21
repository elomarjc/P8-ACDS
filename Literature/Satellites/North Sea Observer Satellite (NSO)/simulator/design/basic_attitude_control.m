% Design of an attitude controller

% State space model of controller
A_att_c = zeros(3,3);
B_att_c = zeros(3,6);
C_att_c = zeros(3,3); 
D_att_c = [0.06555*eye(3) zeros(3)];
