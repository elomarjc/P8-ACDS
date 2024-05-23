clear all
close all
clc
% This file receives the satellite, camera and orbit data, and derives from 
% them the requirements for:
%   - Target tracking : The model with variables definitions can be found in Cubesat ADCS 2023 - Project library
%   - Detumbler
%   - Switching from Nadir to Target tracking.
%   - Nadir pointing


% This file works for:
% - All U satellites satellites

%% Input data: CHANGE THE INPUT DATA HERE, RESPECT THE SELECTED UNITS (International Units)
% Satellite data --------------------------------------------------------------------------------
satellite_type = 2; % 1 = 1U, 2 = 2U
if satellite_type == 1
satellite_inertia_matrix = diag([0.0020,0.0020,0.0020]); % kg*m^2, Satellite moment of inertia, gotten from AAUSAT3
else
    satellite_inertia_matrix = diag([0.0088,0.0088,0.0044]);
end
% Camera information (16 mm Telephoto Lens for RaspberryPi HQ camera) --------------------------------------

%focal_lens  proposed: [16, 25, 50] not decided
camera_focal_length = 16/1000; % m
camera_resolution = [4056,3040];% pixels
camera_pixel_size = 1.55e-6; % m
camera_shutter_speed = 300/1e6; %s, time it takes to take a picture. source: https://forums.raspberrypi.com/viewtopic.php?t=323983


% Orbit information (gotten from aausat 4) Sun Synchoronous orbit
%AAUSAT4_data = [397,502,7.63,5652,98.2,6829.5]
%orbits proposed = [400, 550]*1e3; % m 
orbit_perigee = 500*1e3; % m
orbit_apogee = 600*1e3; % m
orbit_mean_speed = 7.63*1e3; % m/s
orbit_period = 5652; %s
orbit_inclination = 98.2/180*pi;  %rad
orbit_semimajor_axis = 6829.5*1e3; % m

% Element to photograph information: Aalborg
target_name = "Aalborg";
target_area = [60,46]*1e3; % m


% Extra data:
earth_radius = 6371*1e3; %m
earth_average_magnetic_field_surface =  50e-6; % T
percentage_seen = 0.2; % USED FOR SETTLING TIME FROM NADIR TO TARGET TRACKING REQUIREMENT
detumbler_maximum_detumbling_speed = 2; % Hz source: Attitude Determination and Control System for AAUSAT 3
detumbler_end_detumbling_speed = 9.1e-3; % rad/s, (0.3 deg/s in each axis), source: AAUSAT6 Master thesis 2015
detumbler_settling_time = 3; % orbits;
%% Target tracking pass useful variables information 

pass_distance_target_to_satellite = orbit_perigee;% km, average orbiting distance to the surface at this point of the orbit. Equivalent as it assummes a circular orbit.
    % Currently using the worst estimate, to be substituted by the real expected altitude when passing. (defined by orbit)
pass_satellite_radius = earth_radius+pass_distance_target_to_satellite; % m, average orbiting distance to the center of earth at this point of the orbit
pass_earth_angle_start = acos(earth_radius/pass_satellite_radius); % rad, starting angle for target tracking: as soon as it is on view. Angle from the target to the satellite with the origin at the earth's centre.
pass_angular_speed = orbit_mean_speed/(pass_satellite_radius); % rad/s, nadir average angular speed of the satellite with origin at the earth's center (without considering the pointing)
pass_angle_to_nadir_start = pi/2-pass_earth_angle_start; % rad, the starting angle for target tracking from the point of view of the satellite from nadir.

%% Target tracking graphics:
% Assumming t0 = 0, the satellite goes from left to right (2D).
% Nomenclature: 
%  _s : symbolic function as function of time
%  _t : vector of data
%  _f : matlab function handle
syms t

simulation_time_end = 2*pass_earth_angle_start/pass_angular_speed; % s, how long the pass takes.
simulation_time_t = linspace(0,simulation_time_end,1000); % s, vector with time simulation
simulation_Ts = simulation_time_end/1000;
% Satellite earth angle and distance to aalborg
simulation_earth_angle_s =  pass_earth_angle_start - pass_angular_speed*t; % rad, 0 @ target, and positive counterclockwise. Angle from center of the earth, same as pass_earth_angle_start
simulation_distance_target_to_satellite_s =  sqrt(earth_radius^2+pass_satellite_radius^2 -2*earth_radius*pass_satellite_radius*cos(simulation_earth_angle_s)); % m , self defined

% Satellite references: Angle to nadir and its derivatives
simulation_angle_to_nadir_s = asin(sin(simulation_earth_angle_s)./simulation_distance_target_to_satellite_s*earth_radius); % rad, angle to nadir viewed from the satellite
simulation_angular_speed_to_nadir_s =  pass_angular_speed*pass_satellite_radius*cos(simulation_angle_to_nadir_s)./simulation_distance_target_to_satellite_s; % rad/s, derivative of simulation_angle_to_nadir_s
simulation_angular_acceleration_to_nadir_t =  diff(simulation_angular_speed_to_nadir_s); % rad/s^2, second derivative of simulation_angle_to_nadir_s

% Torque and energy
simulation_torque_required_t = simulation_angular_acceleration_to_nadir_t*max(max(satellite_inertia_matrix));
simulation_energy_t = int(simulation_torque_required_t,t,0,t);

% Transform into function handles
simulation_earth_angle_f  = matlabFunction(simulation_earth_angle_s);
simulation_distance_target_to_satellite_f = matlabFunction(simulation_distance_target_to_satellite_s);
simulation_angle_to_nadir_f = matlabFunction(simulation_angle_to_nadir_s);
simulation_angular_speed_to_nadir_f = matlabFunction(simulation_angular_speed_to_nadir_s);
simulation_angular_acceleration_to_nadir_f = matlabFunction(simulation_angular_acceleration_to_nadir_t);
simulation_torque_required_f = matlabFunction(simulation_torque_required_t);
simulation_energy_f = matlabFunction(simulation_energy_t);


% Plotting the target tracking
f1 = figure(1);
clf(f1);
subplot(2,2,1);
hold on
plot(simulation_time_t,simulation_angle_to_nadir_f(simulation_time_t)*180/pi,"Color",[255,165,0]/255,"DisplayName","\theta(t)","LineWidth",2), fontsize(15,"points");
plot([0,simulation_time_end],[pass_angle_to_nadir_start,pass_angle_to_nadir_start]*180/pi,"b--","DisplayName","\theta_0");
xlim([min(simulation_time_t),max(simulation_time_t)]);
xlabel("Time [s]");
ylabel("Angle [deg]");
title("Satellite angle reference from Nadir");
legend
grid on
subplot(2,2,2);
hold on
plot(simulation_time_t,simulation_distance_target_to_satellite_f(simulation_time_t)*1e-3,"g","DisplayName","d(t)","LineWidth",2), fontsize(15,"points");
plot([0,simulation_time_end],[orbit_perigee,orbit_perigee]*1e-3,"b--","DisplayName","perigee");
xlim([min(simulation_time_t),max(simulation_time_t)]);
xlabel("Time [s]");
ylabel("Distance [km]");
title("Satellite to Target Distance");
legend
grid on
subplot(2,2,3);
hold on
plot(simulation_time_t,simulation_angular_speed_to_nadir_f(simulation_time_t),"b","DisplayName","w_\theta(t)","LineWidth",2), fontsize(15,"points");
xlim([min(simulation_time_t),max(simulation_time_t)]);
xlabel("Time [s]");
ylabel("Angular Speed [rad/s]");
title("Satellite angular speed reference from Nadir");
legend
grid on
subplot(2,2,4);
hold on
plot(simulation_time_t,simulation_angular_acceleration_to_nadir_f(simulation_time_t)*1e3,"Color",[	179, 206, 229]/255,"DisplayName","\alpha_\theta(t)","LineWidth",2), fontsize(15,"points");
xlim([min(simulation_time_t),max(simulation_time_t)]);
xlabel("Time [s]");
ylabel("Angular Accel [mrad/s^2]");
title("Satellite angular acceleration reference from Nadir");
legend
grid on

% Torque and energy
figure(2);
clf(2);
subplot(1,2,1);
hold on
grid on
title("Torque applied to statellite");
plot(simulation_time_t,simulation_torque_required_f(simulation_time_t)*1e9,"Color",[	179, 206, 229]/255,"LineWidth",2,"DisplayName","Applied Torque"), fontsize(15,"points");
xlabel("Time [s]");
ylabel("Torque[nNm]");
legend
subplot(1,2,2);
hold on
grid on
title("Energy supplied");
plot(simulation_time_t,simulation_energy_f(simulation_time_t)*1e6,"LineWidth",2,"Color",[200,200,51]/255,"DisplayName","Supplied Energy"), fontsize(15,"points");
xlabel("Time [s]");
ylabel("Energy[uJ]");
legend

%% Worst target tracking requirement:

fprintf("#########################################################\n");
fprintf("#                Target Tracking mode                   #\n");
fprintf("#########################################################\n");
fprintf("Target tracking requirement:\n");

camera_field_of_view_per_pixel = 2*atan(camera_pixel_size/2/camera_focal_length); % rad, field of view of the camera per pixel
%camera_field_of_view = camera_resolution*camera_field_of_view_per_pixel; % rad, field of view of the camera in general
camera_field_of_view = atan(21/500)*ones(2,1);
fprintf("Camera FoV: %.2f by %.2f degrees.\n",camera_field_of_view(1)*180/pi,camera_field_of_view(2)*180/pi); 


pass_area = 2*pass_distance_target_to_satellite*tan(camera_field_of_view/2); % m, the earth surface that can be seen from the satellite (2 dimensions).
fprintf("Minimum camera area: %.2f x %.2f km.\t %s area: %.2f x %.2f km\n",pass_area(1)*1e-3,pass_area(2)*1e-3,target_name,target_area(1)*1e-3,target_area(2)*1e-3);


% Margin of error, considering that aalborg must be contained inside the picture: Longest aalborg dimension vs smallest camera dimension
%pass_slack_area_plus_minus = (min(pass_area)-max(target_area)) / 2; % m, leftover distance to each side of the picture to still have the whole target in frame
pass_slack_area_plus_minus = (21000-7000)/2;
%requirement_angle_accuracy = pass_slack_area_plus_minus/min(pass_area)*min(camera_field_of_view); % rad, the derived angle accuracy required at the top of the path.
requirement_angle_accuracy = pass_slack_area_plus_minus/min(21000)*min(camera_field_of_view);
    % Note that this is a function of the distance between the satellite and the surface (simulation_distance_target_to_satellite_s) so this angle is actually bigger when not on top.
fprintf("Worst scenario pointing requirement (right over %s): %.2f degrees.\n",target_name, requirement_angle_accuracy*180/pi);



%% Torque requirement: What is the maximum torque required for this operation?
pass_angular_speed_to_nadir_maximum = pass_angular_speed*pass_satellite_radius/(pass_satellite_radius-earth_radius); % rad/s, maximum reference angular speed
pass_angular_accel_to_nadir_maximum = max(simulation_angular_acceleration_to_nadir_f(simulation_time_t)); % rad/s^2, maximum angular acceleration reference for target tracking
requirement_torque = pass_angular_accel_to_nadir_maximum*max(max(satellite_inertia_matrix)); % Torque required in the worst case scenario

fprintf("\n----------------------------------------------------------\n");
fprintf("Slew rate requirement:\n");
fprintf("Maximum required slew rate for reference tracking: %.3f mrad/s \n",pass_angular_speed_to_nadir_maximum*1e3);
fprintf("Equivalent torque for reference tracking: %.3f uNm \n",requirement_torque*1e6);

%% Jitter and drift requirement:
% Jitter: The maximum the satellite can jitter while taking a picture, such that it doesn't blur is:
pass_jitter_maximum = camera_field_of_view_per_pixel/camera_shutter_speed; % rad/s, maximum angular speed (w.r.t. the target) the satellite can have before the pixel changes while taking a pic. 
requirement_jitter = pass_jitter_maximum/2/pi; % Hz, maximum angular speed (w.r.t. the target) the satellite can have before the pixel changes while taking a pic. 

fprintf("\n----------------------------------------------------------\n");
fprintf("Jitter and drift requirement:\n");

fprintf("Maximum allowable jitter without blur: %.3f rad/s\n",requirement_jitter);
fprintf("Maximum reference angular speed: %.3f rad/s\n",pass_angular_speed_to_nadir_maximum);
fprintf("Relative magnitude: %.2f dB\n",20*log10(pass_jitter_maximum/pass_angular_speed_to_nadir_maximum));
% Drift: 
fprintf("\nDrift is of no interest in this operation, as it only changes where aalborg is in the picture\n");


%% Settling time requirement:
% The starting time of this operation is defined as the time when aalborg is first in line of sight, so at pass_earth_angle_start angle.
% The ending time shows more slack on its definition: That time when the earth distances appear as 20% of their value by angle effects (not distance to satellite).

settling_angle_to_nadir = asin(percentage_seen); % rads. Angle wrt nadir
[aux,index] = min(simulation_angle_to_nadir_f(simulation_time_t(1:end/4))-settling_angle_to_nadir);
requirement_settling_time = simulation_time_t(index);
settling_time_angle_to_nadir =  simulation_angle_to_nadir_f(requirement_settling_time);
clear aux index

fprintf("\n----------------------------------------------------------\n");
fprintf("Settling time requirement:\n");
fprintf("It takes %.2f s to reach the %.1f %% view.\n",requirement_settling_time,percentage_seen*100);

%% Requirement overshoot from nadir to tracking (for target tracking controller)
% The target must always be in view while shifting from nadir to target tracking. Shift occurs when Aalborg is visible, and finishes
pass_starting_distance = simulation_distance_target_to_satellite_f(0); % m, distance from the satellite to the target at the controller shift
settling_area_seen = 2*pass_starting_distance*tan(camera_field_of_view/2); % m, Theoretical Area in view at the controller shift (not the real one, notice that aalborg is not perpendicular)

% Margin of error, considering worst case with aaborg in picture: Longest aalborg dimension vs smallest camera dimension
settling_slack_area_plus_minus = (min(settling_area_seen)-max(target_area))/2; % m
requirement_overshoot_angle = settling_slack_area_plus_minus/min(settling_area_seen)*min(camera_field_of_view); % rad, overshoot angle for the angle step reference of pass_angle_to_nadir_start


fprintf("\n----------------------------------------------------------\n");
fprintf("nadir to target tracking overshoot requirement (for target tracking controller):\n");
fprintf("Worst scenario overshoot slack: %.2f degrees.\n", requirement_overshoot_angle*180/pi);
fprintf("Nadir angle reference step: %.2f\n",pass_angle_to_nadir_start*180/pi);
fprintf("Maximum percentage overshoot: %.1f %%\n",requirement_overshoot_angle*100/pass_angle_to_nadir_start);

%% Requirements for Detumbler:
% Look at other projects to see how it goes
fprintf("\n\n\n");
fprintf("#########################################################\n");
fprintf("#                      Detumbler mode                   #\n");
fprintf("#########################################################\n");
fprintf("Detumbler maximum recover speed: %.1fHz == %.2f rad/s",detumbler_maximum_detumbling_speed,detumbler_maximum_detumbling_speed*2*pi); 
fprintf("\n----------------------------------------------------------\n");
fprintf("The controller must be able to reduce the tumbling to %.2f mrad/s \n",detumbler_end_detumbling_speed*1e3); 
fprintf("Max angular speed target tracking can recover from: %.2f mrad/s\n",1e3*pass_angular_speed_to_nadir_maximum); % Using the biggest speed reference from the orbit
if detumbler_end_detumbling_speed<pass_angular_speed_to_nadir_maximum
    sentence = "Target tracking can take over after detumbler";
else
    sentence = "Target tracking can NOT take over after detumbler";
end
fprintf("%s",sentence);

% Observed in space: 1.5Hz @ 20 mW consumption, https://studentspace.aau.dk/aausat3/index.php?n=Main.First100Days

fprintf("\n----------------------------------------------------------\n");
fprintf("Settling time requirement: %.1f  orbits\n",detumbler_settling_time); % Same as AAUSAT3


%% Requirements for Nadir pointing:

% Constant slew rate:
fprintf("#########################################################\n");
fprintf("#                          Nadir mode                   #\n");
fprintf("#########################################################\n");

nadir_slew_rate = 2*pi/orbit_period; % rad/s
fprintf("Required slew rate: %.2f mrad/s\n",nadir_slew_rate*1e3);

