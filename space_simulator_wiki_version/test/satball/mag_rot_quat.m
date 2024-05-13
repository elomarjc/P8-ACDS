clc
close all
clear all

% Magnetic field at AAU Student Space
% Coordinates: 57.014361°N 9.9862139°E 
% Magnetic field vector [19-05-2016]: 
% Declination: 2.5911°    [+ East, - West]
% Inclination: 70.8975°   [+ Down, - Up]

% Rotation matrix for adjusting magnetic field to local coordinates
Rx = [1     0               0;          % No rotation around x-axis
      0     cosd(-70.8975)  -sind(-70.8975);  % Rotation around y-axis
      0     sind(-70.8975)  cosd(-70.8975)]; % Rotation around y-axis

% Ry = [cosd(-70.8975) 0 -sind(-70.8975);
%       0              1  0;
%       sind(-70.8975) 0  cosd(-70.8975)];

% Rotation matrix for adjusting magnetic field to local coordinates
Rz = [cosd(-2.5911) -sind(-2.5911) 0;    % Rotation around z-axis
      sind(-2.5911)  cosd(-2.5911) 0;    % Rotation around z-axis
      0              0             1];   % No rotation around z-axis

% Combine the rotation matrices
Rzx = Rz*Rx;

% Convert rotation matrix to quaternion representation
q_zx = rotm2quat(Rzx); % q0 as real part
q_zx_q4 = [q_zx(2) q_zx(3) q_zx(4) q_zx(1)]' % q4 as real part

% Magnetic field vector
% v_m = [0 1 0]*Rzx
v_m = [0.1633 0.1529 0.9746]';

% Normalize the magnetic field vector
v_hat_m = v_m/norm(v_m)

% Rotation matrix for another set of adjustments
Ry = [cosd(180) 0 -sind(180);    % Rotation around y-axis
      0              1  0;       % No rotation around y-axis
      sind(180) 0  cosd(180)];   % Rotation around y-axis

Rz = [cosd(-90) -sind(-90) 0;    % Rotation around z-axis
      sind(-90)  cosd(-90) 0;    % Rotation around z-axis
      0              0             1];   % No rotation around z-axis
  
% Combine the rotation matrices
Rzx = Rz*Rx;

% Convert rotation matrix to quaternion representation
q_zx = rotm2quat(Rzx); % q0 as real part
q_zx_q4 = [q_zx(2) q_zx(3) q_zx(4) q_zx(1)]' % q4 as real part
