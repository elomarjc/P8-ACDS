close all
clc
%This file build upon the requirements_calculator file, run it first.
% This file can take:
% 1U, 1.5U and 2U satellites


% Satellite data:
distance_to_area_centroid = 0.1*satellite_type/2; % So worst case is if all force is at one end of the satellite
gravity_center_deviation = 0.2*(satellite_type==1)+0.02*(satellite_type==1.5)+0.045*(satellite_type==2); % m, +-2cm for a cubesat, with the center of mass at the other end.
    % notice that for the 1.5 U the worst case is still at +-2cm, while for 2U is instead on the z axis (+-4.5cm).
satellite_area = 0.1^2*satellite_type*(satellite_type < 2) + 0.1^2*(satellite_type==2); % m^2, notice that, as a 2U, the worst case is around the z axis.

drag_coefficient_cube = 1.05; % https://study.com/academy/lesson/drag-coefficient-overview-equation.html
%drag_coefficient_cube = 2.5; % Worst case https://ntrs.nasa.gov/api/citations/20110007876/downloads/20110007876.pdf

satellite_magnetic_moment =  0.0091*satellite_type; % Am^2. Source: NASA paper / wertz book. it suggest 1Am^2 for a 110W satellite, so 0.0091 Am^2/W

    
fprintf("#########################################################\n");
fprintf("#                Disturbance estimation                 #\n");
fprintf("#########################################################\n");


%% Solar radiation
fprintf("\n----------------------------------------------------------\n");
fprintf("Solar radiation:\n");
% Values: https://aaudk.sharepoint.com/sites/ES8-Masters/Delte%20dokumenter/General/Literature/NASA%20material%20for%20CubeSat%20development/ADCS%20design%20guide%20-%20NASA.pdf?CT=1710936171463&OR=ItemsView
solar_constant = 1367; %W/m^2, average solar radiation from earth distance
speed_of_light = 3e8; % m/s
reflection_factor = 1; % e[0,1] Worst scenario
solar_radiation_torque = @(ang_incidence,angle_cA_to_cG) abs((distance_to_area_centroid+gravity_center_deviation))*satellite_area*(1+reflection_factor)*(solar_constant/speed_of_light)*cosd(ang_incidence)*sind(angle_cA_to_cG);
% Angles in degrees
% Maximum solar radiation torque: The cm is at the limit of the specification, the light incedes perpendiculary 
max_solar_torque = solar_radiation_torque(0,90);
fprintf("Worst solar radiation pressure: %.2f nNm\n",max_solar_torque*1e9);
% Value consistent with literature: https://elib.dlr.de/128377/1/Disturbance_Torques_Upon_GRACE_MichaelTolstoj.pdf
% When does this worst happen? Depends on the orbit of the satellite
%% Aerodynamic drag
fprintf("\n----------------------------------------------------------\n");
fprintf("Atmospheric drag:\n");
% Modeling disturbances influencing an Earth-orbiting satellite

density_atm_low= 1.585e-12; %kg/m3 @ 450
density_atm_high = 6.967e-13;%kg/m3 @ 500
%Source: https://aaudk.sharepoint.com/sites/ES8-Masters/Delte%20dokumenter/General/Literature/Satellites/AAUSAT6/2023%20-%20BP%20-%20CubeSat%20ADCS.pdf?CT=1710937185865&OR=ItemsView

aerodynamic_torque = @(density_low,density_high)  abs(distance_to_area_centroid+gravity_center_deviation)*1/2*abs(density_low-density_high)*drag_coefficient_cube*(satellite_area/2)*(orbit_mean_speed)^2;

max_aerodynamic_torque = aerodynamic_torque(density_atm_low,density_atm_high);
fprintf("Worst aerodynamic torque: %.2f nNm\n",max_aerodynamic_torque*1e9);

%% magnetic torque
fprintf("\n----------------------------------------------------------\n");
fprintf("Magnetic drag:\n");
%earths_max_magnetic_field = 60e-6; % T source: Wertz, fig 5.1 at around 400km. Remember to multiply by 2, as at the pole is x2 the equator(graph)
magnetic_field_at_orbit = earth_average_magnetic_field_surface*(earth_radius/pass_satellite_radius)^3 % T
max_magnetic_torque = satellite_magnetic_moment*magnetic_field_at_orbit;
fprintf("Worst magnetic torque: %.2f nNm\n",max_magnetic_torque*1e9);
fprintf("It is presummed that it can be predicted and thus removed up to 90%% by magnetorquers\nEquivalent torque: %.2f nNm\n",0.1*max_magnetic_torque*1e9);

%% Gravity torque
fprintf("\n----------------------------------------------------------\n");
fprintf("Gravity effect:\n");
earth_mass = 5.972e24; %kg
earth_radii = 6.371e6; % m
gravitational_constant = 6.67e-11; % N m^2/kg^2
max_torque_gravity = 3*earth_mass*gravitational_constant/2/(orbit_perigee+earth_radii)^3*norm(satellite_inertia_matrix);
fprintf("Worst gravity torque: %.2f nNm\n",max_torque_gravity*1e9);


%% Overall result
fprintf("\n----------------------------------------------------------\n");
fprintf("Result:\n");
worst_torque = sqrt(max_solar_torque^2+max_aerodynamic_torque^2+(max_magnetic_torque)^2+max_torque_gravity^2);
worst_torque_dipole = sqrt(max_solar_torque^2+max_aerodynamic_torque^2+(max_magnetic_torque/10)^2+max_torque_gravity^2);

%Worst case scenario:
fprintf("Worst case scenario: %.2f nNm\n",worst_torque*1e9);
fprintf("With dipole moment reduction of 90%%: %.2f nNm\n",worst_torque_dipole*1e9);
fprintf("Maximum required torque for reference tracking: %.2f nNm \n",requirement_torque*1e9);
