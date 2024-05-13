%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is used to find the worst case gravity
% gradient torque. It is assumed that the inertia
% matrix is as found the report. From the combined
% inertia matrix of the satellite body and solar
% arrays the direction yielding the larges gravity
% gradient torque is found.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R_e      = 6378000; % Radius of the Earth [m]
orbit    = 500000;  % Orbit Height [m]
my       = 3.986e14;% Gravitational constant
m_flap   = 0.5;     % mass of solar array [kg]
m_sat    = 2.5;     % mass of satellite body [kg]

GOM      = [0.05,0.05,0.15]'; % placement of geometric center in SBRF (x,y,z) [m]
COM_flap = [0.05, 0.05, 0]';  % placement of solar array CoM in SBRF (x,y,z) [m]
COM_sat  = GOM;               % placement of satellite body CoM in SBRF (x,y,z) [m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Used when worst worst case is calculated
%  - one mass at each end of the satellite
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ r_offset = [0 0 0.02]';     % offset of satellite body CoM in SBRF (x,y,z) [m]
% $$$ COM_sat = COM_sat+r_offset;
% $$$ x1 = [0.05 0.05 30]';       % position of mass 1 in SBRF (x,y,z) [m]
% $$$ x2 = [0.05 0.05 0]';        % position of mass 2 in SBRF (x,y,z) [m]
% $$$ 
% $$$ m = pinv([x1 x2])*COM_sat*2.5; % calculates masses
% $$$ 
% $$$ %calculation of Inertia Matrix
% $$$ I_xx = (x1(2)^2+x1(3)^2)*m(1)+(x2(2)^2+x2(3)^2)*m(2);
% $$$ I_yy = (x1(1)^2+x1(3)^2)*m(1)+(x2(1)^2+x2(3)^2)*m(2);
% $$$ I_zz = (x1(1)^2+x1(2)^2)*m(1)+(x2(1)^2+x2(2)^2)*m(2);
% $$$ I_xy = x1(1)*x1(2)*m(1)+x2(1)*x2(2)*m(2);
% $$$ I_xz = x1(1)*x1(3)*m(1)+x2(1)*x2(3)*m(2);
% $$$ I_yz = x1(2)*x1(3)*m(1)+x2(2)*x2(3)*m(2);
% $$$ I_sat = [I_xx -I_xy -I_xz;-I_xy I_yy -I_yz;-I_xz -I_yz I_zz];

%I_sat = diag([0.0208,0.0208,0.0042]);    % Inertia Matrix of satellite body
%I_flap = diag([0.01125 0.01125 0.0225]); % Inertia Matrix of solar array
%I_sat_folded = diag([0.025 0.025 0.005]);% Inertia Matrix of satellite in
                                         % folded state

I_deployed = diag([0.0422 0.0422 0.0283]);
I_folded = diag([0.025 0.025 0.005]);
%COM = (COM_sat*m_sat+COM_flap*m_flap)/(m_sat+m_flap); % combined CoM
COM_deployed = [0.05 0.05 0.125]';
COM_folded = [0.05 0.05 0.15]';

% $$$ % calculates offset of CoM from combined CoM
% $$$ offset_flap = COM-COM_flap;
% $$$ offset_sat = COM-COM_sat;
% $$$ 
% $$$ % calculates Inertia Matrix for combined CoM
% $$$ I_sat_COM = I_sat+m_sat*(offset_sat'*offset_sat*eye(3)-offset_sat*offset_sat');
% $$$ I_flap_COM = I_flap+m_flap*(offset_flap'*offset_flap*eye(3)-offset_flap*offset_flap');
% $$$ I_total = I_flap_COM+I_sat_COM;

% Finds the eigen vectors and eigen values, i.e., the 
% rotation matrix and Inertia matrix in principal axis
[C,I_pri] = eig(I_deployed);
[C_f,I_f] = eig(I_folded);

% This is the earth direction yielding the largest gravity
% gradient. This value assumes abs(I_xx-I_zz) is larger than
% abs(I_xx-I_yy) and abs(I_yy-I_zz)
att_w = [1 0 1]';%[1/2 0 sqrt(1-0.5^2)]';
% rotates earth direction into SBRF 
att = C*att_w;
att = att/norm(att);

att_f = C_f*att_w;
att_f = att_f/norm(att_f);

% calculates the torque
torque_gg = 3*my/(R_e+orbit)^3*cross(att,I_deployed*att);
torque_gg_f = 3*my/(R_e+orbit)^3*cross(att_f,I_folded*att_f);

disp('Max Gravity Gradient Torque - deployed')
norm(torque_gg)

disp('Max Gravity Gradient Torque - folded')
norm(torque_gg_f)

