%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file calculates the inertia matrix of the satellite in two
% different cases.
%
% Case 1: The solar arrays have been deployed.
%         In this case the following approximations have been made:
%         - The satellite body is approximated as a solid brick with mass 
%         M1 = 2.5 kg and the dimensions 0.1x0.1x0.3 m.
%         - The solar arrays are approximated as four solid bricks with 
%         mass M2 = 0.125 kg each were the hight is 2 mm.
%
% Case 2: The solar arrays have not been deployed.
%         - The satellite is approximated as a solid brick with mass M = 3
%         kg and the dimensions 0.1x0.1x0.3 m.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc;

% Constants
M1 = 2.5;
M2 = 0.125;
M = M1+4*M2;
R = 0.30;
% Dimensions of satellite
a = 0.1;
b = 0.1;
c = 0.3;
% Dimensions of solar arrays
a1 = 0.1;
b1 = 0.3;
c1 = 2e-3;

%%%%%%%%%%
% Case 1 %
%%%%%%%%%%

% The inertia matrix of the solid brick around the Center of Mass (CoM) can
% be calculated in the following way:

A = (1/12)*M1*(b^2+c^2);
B = (1/12)*M1*(a^2+c^2);
C = (1/12)*M1*(a^2+b^2);

I1 = [A 0 0;0 B 0;0 0 C];

% The inertia matrix of each of the solar arrays are equal in pairs so two
% inertia matrices are needed here these are:

D = (1/12)*M2*(b1^2+c1^2);
E = (1/12)*M2*(a1^2+c1^2);
F = (1/12)*M2*(a1^2+b1^2);

I2 = [D 0 0;0 E 0;0 0 F];
I3 = [E 0 0;0 D 0;0 0 F];

% In order to find the combined CoM vectors from the origin of the SBRF to
% the CoM of each object are needed these are defined as follows:

p1 = [0.05 0.05 0.15]';
p2 = [0.05 -0.15 0]';
p3 = [0.05 0.25 0]';
p4 = [0.25 0.05 0]';
p5 = [-0.15 0.05 0]';

% The combined CoM can then be calculated as

v_com = (M1*p1 + M2*p2 + M2*p3 + M2*p4 + M2*p5)/(M1 + 4*M2);

% The inertia matrices must be move to the combined CoM. This is done with
% the following equation I' = I + M*(p'*p*eye(3) - p*p').

I_com1 = I1 + M1*((v_com-p1)'*(v_com-p1)*eye(3) - (v_com-p1)*(v_com-p1)');
I_com2 = I2 + M2*((v_com-p2)'*(v_com-p2)*eye(3) - (v_com-p2)*(v_com-p2)');
I_com3 = I3 + M2*((v_com-p3)'*(v_com-p3)*eye(3) - (v_com-p3)*(v_com-p3)');
I_com4 = I2 + M2*((v_com-p4)'*(v_com-p4)*eye(3) - (v_com-p4)*(v_com-p4)');
I_com5 = I3 + M2*((v_com-p5)'*(v_com-p5)*eye(3) - (v_com-p5)*(v_com-p5)');

% The inertia matrix of the satellite around the combined CoM is then
disp('Inertia matrix around the conbined CoM, when solar arrays are deployed')
I = I_com1 + I_com2 + I_com3 + I_com4 + I_com5

%%%%%%%%%%
% Case 2 %
%%%%%%%%%%

% The inertia matrix of the solid brick around the Center of Mass (CoM) can
% be calculated in the following way:

F = (1/12)*M*(b^2+c^2);
G = (1/12)*M*(a^2+c^2);
H = (1/12)*M*(a^2+b^2);

disp('Inertia matrix around CoM when solar arrays are not deployed')
I = [F 0 0;0 G 0;0 0 H]
