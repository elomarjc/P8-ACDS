%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file calculates the worst case aerodynamic
% torque. It is assumed that the satellite consists
% of two bodies - the solar array and satellite
% body. The solar array has a fixed CoM where the
% satellite body's CoM is changed within 2cm of the
% geometrical center. First the combined CoM for
% the satellite body and solar array is calculated.
% Then the area perpendicular to the velocity
% vector is calculated and finally the aerodynamic
% force and torque are calculated.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
printEPS = true;              % print the plot to a file if true
Cd       = 2;                 % Drag coefficient
rho      = 6.0e-12+7.3e-13;             % density of the atmosphere
v        = 7613;              % air speed [m/s]
m_sat    = 2.5;               % mass of satellite body [kg]
m_flap   = 0.5;               % mass of solar array [kg]
GOM      = [0.05,0.05,0.15]'; % position of the Geometric center in SBRF (x,y,z) [m]
COM_flap = [0.05,0.05,0]';    % position of solar array CoM in SBRF (x,y,z) [m]
z_max    = 100;                % maximum value of the z-axis component of
                              % the velocity vector
R_ref    = [1,1,-z_max];      % air velocity unit vector start direction in SBRF (x,y,z) [m]


% dimensions of satellite body
x = [0.10 0 0]';              % dimension of satellite body x-axis in SBRF (x,y,z) [m]
y = [0 0.10 0]';              % dimension of satellite body y-axis in SBRF (x,y,z) [m]
z = [0 0 0.30]';              % dimension of satellite body z-axis in SBRF (x,y,z) [m]

% dimensions of solar arrays
f1_x = x;                     % dimension of solar array 1 x-axis in SBRF (x,y,z) [m]
f1_y = [0 -0.30 0]';          % dimension of solar array 1 y-axis in SBRF (x,y,z) [m]
f2_x = [-0.30 0 0]';          % dimension of solar array 2 x-axis in SBRF (x,y,z) [m]
f2_y = y;                     % dimension of solar array 2 y-axis in SBRF (x,y,z) [m]

% arrays used to store max torques at each CoM
torque_a_all = [];

% for loop used to change position of CoM with a constrained lenght of 2cm
% arrays used to store torques at each air velocity direction
torque_a_no_flap = [];
torque_a = [];

angle_a = [];

% calculates CoM of the satellite body and combined CoM
COM_sat = GOM;
COM = (COM_sat*m_sat+COM_flap*m_flap)/(m_sat+m_flap);

% constants used in the projection matrix
A = R_ref(1);
B = R_ref(2);
% for loop used to change air velocity direction
for C = R_ref(3):0.1:z_max;
  % calculates the air velocity unit vector
  v_hat = [A,B,C]';
  v_hat = v_hat/norm(v_hat);
  
  % calculates the projection matrix
  P = 1/(A^2 + B^2 + C^2)*[B^2+C^2, -A*B, -A*C;-A*B,A^2+C^2,-B*C;-A*C,-B*C A^2+B^2];
  
  % projects the satellite body and solar array dimension onto the
  % plane ortogonal to the air velocity vector
  xp = P*x;
  yp = P*y;
  zp = P*z;
  f1_px = P*f1_x;
  f1_py = P*f1_y;
  f2_px = P*f2_x;
  f2_py = P*f2_y;
  
  % projects the the z-axis of the satellite body in the plane onto the
  % solar array axes in the plane 
  f1_pxp = dot(zp,f1_px)/norm(f1_px)^2*f1_px;
  f1_pyp= dot(zp,f1_py)/norm(f1_py)^2*f1_py;
  f2_pxp = dot(zp,f2_px)/norm(f2_px)^2*f2_px;
  f2_pyp = dot(zp,f2_py)/norm(f2_py)^2*f2_py;
  
  % determines if the shadow of the satellite body z-axis is greater than the solar array
  if(norm(f1_pyp)> norm(f1_py))
    f1_pyp = f1_py;
  end
  if(norm(f2_pxp)>norm(f2_px))
    hep = f2_px;
  else
    hep = f2_pxp;
  end
  
  % calculates the shadow area in the plane of the different solar arrays. The shadow area
  % is either a triangle or a square - a triangle.
  if(norm(f1_pxp) < norm(f1_px))
    A_shad_1 = norm(cross(xp,f1_pyp))-0.5*norm(cross(f1_pxp,f1_pyp));
  else
    A_shad_1 = 0.5*norm(cross(xp,f1_pyp));
  end
  if(norm(f2_pyp) < norm(f2_py))
    A_shad_2 = norm(cross(yp,hep))-0.5*norm(cross(f2_pyp,hep));
  else
    A_shad_2 = 0.5*norm(cross(yp,hep));
  end
  
  % calculates the combined shadow area in the plane
  A_shadow = A_shad_1+A_shad_2;
  % calculates the area of the four solar arrays in the plane
  A_flap = norm(cross(f1_px,f1_py))*4;
  % calculates the projected area of the satellite body.
  A_sat = norm(cross(xp,yp))+norm(cross(yp,zp))+norm(cross(xp,zp));
  
  A_tot = A_sat+A_flap-A_shadow;
  % calculates the aerodynamic force on the satellite where the solar
  % array is deployed
  F_a = -0.5*rho*Cd*v^2*A_tot;
  % calculates the aerodynamic force on the satellite where the solar
  % array is folded
  F_a_no_flap = -0.5*rho*Cd*v^2*A_sat;
  % torques with folded solar array at different attitudes are "stored"
  torque_a_no_flap = [torque_a_no_flap,norm(cross(GOM-COM_sat,F_a_no_flap*v_hat))];
  % torques with deployed solar array at different attitudes are "stored"
  torque_a = [torque_a,norm(cross(GOM-COM,F_a*v_hat))];
  if(C<0)
    angle_a = [angle_a,-acos(dot([1 1 0]'/norm([1 1 0]'),v_hat))];
  else
    angle_a = [angle_a,acos(dot([1 1 0]'/norm([1 1 0]'),v_hat))];
  end
end

figure()
plot([-90 angle_a.*(180/pi) 90],[0 torque_a 0])
grid on
xlabel('Angle between solar array plane and velocity vector [\circ]')
ylabel('Torque due to aerodynamic drag [Nm]')
title('Aerodynamic Torque')

if(printEPS)
    print -depsc aerodynamic_torque_vs_attitude.eps
end