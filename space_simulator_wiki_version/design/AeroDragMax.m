
if exist('RunFile','var')==0
    clc
    close all
    clear
    WCDTorque
end

%% Init
if exist('RunFile','var')==1

R_ref    = [1,0,1];           % air velocity unit vector start direction in SBRF (x,y,z) [m]

% Satellite dimensions
x = [SatWidth 0 0]';              % dimension of satellite body x-axis in SBRF (x,y,z) [m]
y = [0 SatLength 0]';              % dimension of satellite body y-axis in SBRF (x,y,z) [m]
z = [0 0 SatHeigth]';              % dimension of satellite body z-axis in SBRF (x,y,z) [m]

% arrays used to store max torques at each CoM
m_torque = [];            

% arrays used to store torques at each air velocity direction
torque = [];
A2 = [];

% constants used in the projection matrix
A = R_ref(1);
C = R_ref(3);
% for loop used to change air velocity direction
  for B = R_ref(2):0.1:10;
    % calculates the air velocity unit vector
    v_hat = [A,B,C]';
    v_hat = v_hat/norm(v_hat);
    
    % calculates the projection matrix
    P = 1/(A^2 + B^2 + C^2)*[B^2+C^2, -A*B, -A*C;-A*B,A^2+C^2,-B*C;-A*C,-B*C A^2+B^2];
    
    % projects the satellite onto the plane ortogonal to the air velocity vector
    xp = P*x;
    yp = P*y;
    zp = P*z;
    
    A_sat = norm(cross(xp,yp))+norm(cross(yp,zp))+norm(cross(xp,zp))
    F = -0.5*rho*CD*V^2*A_sat;
    % torques at different attitudes are "stored"
    torque = [torque,norm(cross(GoM-CoM,F*v_hat))];
    A2=[A2,A_sat];
  end
  
% maximum torque is found
m_torque = max(torque);


disp('Max Aerodynamic Torque')
max(m_torque)

end
