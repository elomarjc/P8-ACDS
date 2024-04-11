close all
clc
%This file build upon the requirements_calculator file and the disturbances_calculator file, run it first

%% Flywheel sizing:

margin_factor = 1.2; % 20% margin factor



fprintf("#########################################################\n");
fprintf("#                   Flywheel sizing                     #\n");
fprintf("#########################################################\n");

% Maximum torque to apply: ------------------------------

fprintf("\n----------------------------------------------------------\n");
fprintf("Maximum torque required:\n");
% Max ref torque: requirement_torque
% Max dist torque: worst_torque

fprintf("Maximum reference torque: %.2f nNm\n",requirement_torque*1e9);
fprintf("Maximum disturbance torque: %.2f nNm\n",worst_torque*1e9);

fprintf("Nadir to target pointing: Less than target pointing:\nIt has to correct %.2f degrees in %.2fs, which is less than the %.2f over the target from 250 to 400s from the pointing .\n",pass_angle_to_nadir_start*180/pi,requirement_settling_time,(simulation_angle_to_nadir_f(250)-simulation_angle_to_nadir_f(400))*180/pi); 

flywheel_torque_required = requirement_torque+worst_torque*1.2; % The total required torque is the tracking maximum torque plus the disturbances applied
fprintf("\nFlywheel torque authority for disturbance rejection + tracking: %.2f nNm\n",flywheel_torque_required*1e9);
% Maximum energy storage required
fprintf("\n----------------------------------------------------------\n");
fprintf("Momentum storage required:\n");
point_tracking_max_energy = max(simulation_energy_f(simulation_time_t)); % Found from integration of the torque in the curves
fprintf("Point tracking maximum energy storage: %.2f uJ\n",point_tracking_max_energy*1e6);
nadir_pointing_max_energy = orbit_period*worst_torque/4; % Integration over one gorth orbit of the disturbances, as they are mainly cyclic
nadir_pointing_max_energy_dipole = orbit_period*worst_torque_dipole/4; % Integration over one gorth orbit of the disturbances, as they are mainly cyclic
fprintf("Nadir pointing maximum energy storage: %2.f uJ. If magnetic dipole corrected: %2.f uJ\n",nadir_pointing_max_energy*1e6,nadir_pointing_max_energy_dipole*1e6);
nadir_to_pointing_max_energy = 1/2*max(max(satellite_inertia_matrix))*(simulation_angular_speed_to_nadir_f(requirement_settling_time))^2 -1/2*max(max(satellite_inertia_matrix))*nadir_slew_rate^2; % Energy difference between nadir pointing and the settled target tracking
fprintf("Nadir to target pointing energy difference: %.2f uJ\n",nadir_to_pointing_max_energy*1e6);
flywheel_max_energy = point_tracking_max_energy+nadir_pointing_max_energy+nadir_to_pointing_max_energy;
flywheel_max_energy_dipole = point_tracking_max_energy+nadir_pointing_max_energy_dipole+nadir_to_pointing_max_energy;
fprintf("\nMaximum energy storage needed, then: %.2f uJ, dipole corrected: %.2f uJ\n",flywheel_max_energy*1e6,flywheel_max_energy_dipole*1e6);


%% Required flywheel
material_name = "brass";
material_density =8.73*1e-3*1e6; % kg/m^3
max_flywheel_speed = 1000*pi/30; %rad/s max desired flywheel speed
radius_to_lenght_ratio = 3; % relation between height and radius

flyhweel_required_inertia = flywheel_max_energy/max_flywheel_speed^2*2 % kg*m^2
flyhweel_required_inertia =     0.38E-6;
flywheel_required_radius = (radius_to_lenght_ratio*flyhweel_required_inertia*2/pi/material_density)^(1/5); % m
flywheel_required_height = flywheel_required_radius/radius_to_lenght_ratio; % m
flywheel_required_mass = flywheel_required_height*pi*material_density*flywheel_required_radius^2; % kg
fprintf("Using a radius to height ratio of : %.1f\n",radius_to_lenght_ratio);
fprintf("Required flywheel mass: %.2f g\n",flywheel_required_mass*1e3);
fprintf("Radius: %.2f mm \t Length: %.2f mm\n",flywheel_required_radius*1e3,flywheel_required_height*1e3);


%% Flywheel simulation with given data:
%motor_R = 4.44; %ohms
%motor_L = 0.12e-3; %H
%motor_Kt = 1.81e-3;% N*m/A
%motor_Ke = 1.81e-3;% N*m/A
%motor_b = 63.90*1e-9; %Nms
%motor_I = 0.381e-6; % kg m2;
%motor_t_mech = 32.1e-3; %s
%s = tf('s');
%G_w_V = motor_Kt/((motor_R+motor_L*s)*(motor_I*s+motor_b)+motor_Ke*motor_Kt)
%G_tau_w = 1/(motor_I*s+motor_b);
%Ks = (0.005619*s+0.01347)/s;
%T = feedback(G_w_V*Ks,1);
 

%t = linspace(0,simulation_time_end,10000);
%u = ones(size(t))*worst_torque;
%lsim(G_tau_w,u,t);
%% Employed flywheel previously:
clear maxon
maxon.R = 12.4; %ohms
maxon.L =0.281E-3; %H
maxon.Kt = 2.79E-3; %Nm/A
maxon.Ke = 1/(3450*pi/30); % V*s/rad
maxon.Ishaft = 0.237E-7; % kgm^2
maxon.tmech = 38.2E-3; % s
maxon.b = maxon.Ishaft/maxon.tmech;
tau_step = 1000e-9; %nN
F = tau_step/maxon.Kt;


s = tf('s');
G = maxon.Kt/(maxon.R+s*maxon.L); % A/V
C = (s+100)/s;
L = G*C;
figure
margin(L)
figure(1)
subplot(2,1,1);
T = feedback(L,1);
step(T);
subplot(2,1,2);
step((1-T)*C*tau_step/maxon.Kt)