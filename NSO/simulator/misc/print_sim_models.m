% Script intended for creating eps-files of the Simulink implementations of
% the different models
% info about open_system:
% http://www.mathworks.com/access/helpdesk/help/toolbox/simulink/slref/open_system.html

clc

% All files are printed to outputDir
outputDir = 'sim_model_figs';
systemStr = computer;
%Make dir
answ = 'n';
fprintf('\n\n')
answ	= input('Caveat Emptor! Simulink Model Printer: \n   Use will delete the old output dir!! \n   Confirm action (y/n) [n]  ','s');
if strcmpi(answ,'y')
    fprintf('\n')
	fprintf('\n\t--= Deleting old output dir =--\n')
else
	fprintf('\n')
	fprintf('\n\t--= PROGRAM TERMINATED =--\n')
	return
end
%make fig dir
current_dir = pwd;
if(strcmp(systemStr,'PCWIN'))
  [status, result] = system(strcat('rmdir /S /Q  "',current_dir,'\',outputDir,'"'));
  if status == 1
    fprintf('\t\tNo old fig directory exists')
  end
else
  [status, result] = system(strcat('rm -rf "',current_dir,'/',outputDir,'"'));
  if status == 1
    fprintf('\t\tNo old fig directory exists')
  end
end
fprintf('\n\t--= Making new output dir =--\n')
[status, result] = system(strcat('mkdir "',current_dir,'/',outputDir,'"'));
if status == 1
    fprintf('\t\tCould not create output dir - returning:')
    fprintf(strcat('\n     ',result))
    fprintf('\n\t\tCurrentdir:\n')
    fprintf(strcat('\n\t\t',current_dir))
    return
end
fprintf('\n     --= Printing systems =--\n')
%%%%%%%%%%%%%%%%%%%%%% Print lib blocks %%%%%%%%%%%%%%%%%%%%
% open top lib
topLibPath = '../lib/';
topLib = 'components_lib';
open_system(strcat(topLibPath,topLib))

%%%%%%%%%% Disturbances %%%%%%%%%%
tempLibObj = 'Disturbances/Radiation (I)';
outputName = 'disturbance_radiation.eps';
print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)

tempLibObj = 'Disturbances/Atmospheric';
outputName = 'disturbance_atmospheric.eps';
print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)

tempLibObj = 'Disturbances/Gravity/Gravity Gradient Torque';
outputName = 'disturbance_gravity.eps';
print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)

% $$$ tempLibObj = 'Disturbances/Magnetic residual disturbance';
% $$$ outputName = 'disturbance_magnetic.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)
% $$$ 
% $$$ tempLibObj = 'Disturbances/Exposed area';
% $$$ outputName = 'disturbance_exposed_area.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)
% $$$ 
tempLibObj = 'Disturbances/Disturbances';
outputName = 'disturbance_total.eps';
print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)

%%%%%%%%%% Orbit and Ephemeris %%%%%%%%%%
% $$$ tempLibObj = 'Orbit and Ephemeris/Albedo emulation model';
% $$$ outputName = 'orbit_ephemeris_albedo_emulation.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)
% $$$ 
% $$$ tempLibObj = 'Orbit and Ephemeris/Ephemeris';
% $$$ outputName = 'orbit_ephemeris_ephemeris.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)
% $$$ 
% $$$ tempLibObj = 'Orbit and Ephemeris/Orbit and Magnetic field models';
% $$$ outputName = 'orbit_ephemeris_total.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)
% $$$ 
% $$$ tempLibObj = 'Orbit and Ephemeris/Eclipse model';
% $$$ outputName = 'orbit_ephemeris_eclipse_model.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)
% $$$ 
% $$$ tempLibObj = 'Orbit and Ephemeris/Albedo';
% $$$ outputName = 'orbit_ephemeris_albedo.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)
% $$$ 
% $$$ tempLibObj = 'Orbit and Ephemeris/Magnetic field model';
% $$$ outputName = 'orbit_ephemeris_magnetic_field.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)
% $$$ 
% $$$ tempLibObj = 'Orbit and Ephemeris/(IGRF) Magnetic field';
% $$$ outputName = 'orbit_ephemeris_IGRF.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)
% $$$ 
%%%%%%%%%% Sensor Emulation %%%%%%%%%%
% $$$ tempLibObj = 'Sensor Emulation/Sensor Emulation';
% $$$ outputName = 'sensor_emulation_total.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)
% $$$ 
% $$$ tempLibObj = 'Sensor Emulation/Emulation of temperature sensors';
% $$$ outputName = 'sensor_emulation_temperature_sensor.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)
% $$$ 
% $$$ tempLibObj = 'Sensor Emulation/Emulation of gyros';
% $$$ outputName = 'sensor_emulation_gyros.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)
% $$$ 
% $$$ tempLibObj = 'Sensor Emulation/Sun sensors';
% $$$ outputName = 'sensor_emulation_sun_sensor.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)

%%%%%%%%%% Dynamics and Kinematics
% $$$ tempLibObj = 'Dynamics and Kinematics/Spacecraft Dynamics';
% $$$ outputName = 'sensor_emulation_dynamics_kinematics.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)

%%%%%%%%%% ADS - on-board models%%%%%%%%%%
% $$$ tempLibObj = 'ADS/On-board model/ADS on-board models';
% $$$ outputName = 'on-board_models_on-board_models.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)
% $$$ 
% $$$ tempLibObj = 'ADS/On-board model/sensor model';
% $$$ outputName = 'on-board_models_sentor_model.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)
% $$$ 
% $$$ tempLibObj = 'ADS/On-board model/Albedo model onboard';
% $$$ outputName = 'on-board_models_albedo_model.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)
% $$$ 
% $$$ tempLibObj = 'ADS/On-board model/Ephemeris on-board';
% $$$ outputName = 'on-board_models_ephemeris.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)
% $$$ 
% $$$ tempLibObj = 'ADS/On-board model/Orbit and Magnetic field models';
% $$$ outputName = 'on-board_models_orbit_mag.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)
% $$$ 
% $$$ tempLibObj = 'ADS/On-board model/Magnetometer';
% $$$ outputName = 'on-board_models_magnetometer.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)
% $$$ 
% $$$ tempLibObj = 'ADS/On-board model/Temperature: Convert voltage from ADC into measured temperature';
% $$$ outputName = 'on-board_models_temperature.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)
% $$$ 
% $$$ tempLibObj = 'ADS/On-board model/Gyros';
% $$$ outputName = 'on-board_models_gyro.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)
% $$$ 
% $$$ tempLibObj = 'ADS/On-board model/Albedo projection';
% $$$ outputName = 'on-board_models_albedo_project.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)
% $$$ 
% $$$ tempLibObj = 'ADS/On-board model/Albedo emulation model';
% $$$ outputName = 'on-board_models_albedo_emu.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)
% $$$ 
% $$$ tempLibObj = 'ADS/On-board model/Eclipse model';
% $$$ outputName = 'on-board_models_eclipse.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)
% $$$ 
% $$$ tempLibObj = 'ADS/On-board model/Sun sensors';
% $$$ outputName = 'on-board_models_sun_sensor.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)
% $$$ 
% $$$ tempLibObj = 'ADS/On-board model/Magnetic field model';
% $$$ outputName = 'on-board_models_mag_field.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)

%%%%%%%%%% ADS - Attittude Determination
% $$$ tempLibObj = 'ADS/Attitude Determination/q-method';
% $$$ outputName = 'attitude_determination_q_method.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)
% $$$ 
% $$$ tempLibObj = 'ADS/Attitude Determination/Propagator from Kalman filter';
% $$$ outputName = 'attitude_determination_propagator_kalman.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)
% $$$ 
% $$$ tempLibObj = 'ADS/Attitude Determination/EKF using sun and mag vectors .. and removing bias';
% $$$ outputName = 'attitude_determination_EKF_sun_mag_bias.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)
% $$$ 
% $$$ tempLibObj = 'ADS/Attitude Determination/otoquest';
% $$$ outputName = 'attitude_determination_otoquest.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)
% $$$ 
% $$$ tempLibObj = 'ADS/Attitude Determination/EKF using quaternions';
% $$$ outputName = 'attitude_determination_EKF_quaternion.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)
% $$$ 
% $$$ tempLibObj = 'ADS/Attitude Determination/otoquest with propagator';
% $$$ outputName = 'attitude_determination_otoquest_propagator.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)
% $$$ 
% $$$ tempLibObj = 'ADS/Attitude Determination/EKF using sun and mag vectors';
% $$$ outputName = 'attitude_determination_EKF_sun_mag.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)

%%%%%%%%%% ACS - Controllers %%%%%%%%%%
% $$$ tempLibObj = 'ACS/Controllers/Detumble controller';
% $$$ outputName = 'controller_detumble.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)
% $$$ 
% $$$ tempLibObj = 'ACS/Controllers/LQ-Regulator';
% $$$ outputName = 'controller_LQ_controller.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)

%%%%%%%%%% Actuator Emulation %%%%%%%%%%
tempLibObj = 'Actuator Emulation/Magnetorquer without L';
outputName = 'actuator_emulation_magnetorquer.eps';
print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)

tempLibObj = 'Actuator Emulation/Momentum wheel';
outputName = 'actuator_emulation_momentum_wheel.eps';
print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)

% $$$ tempLibObj = 'Actuator Emulation/PI-controller';
% $$$ outputName = 'actuator_emulation_pi_controller.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)
% $$$ 
% $$$ tempLibObj = 'Actuator Emulation/NSO magnetorquers';
% $$$ outputName = 'actuator_emulation_nso_magnetorquers.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)
% $$$ 
% $$$ tempLibObj = 'Actuator Emulation/NSO mw HW';
% $$$ outputName = 'actuator_emulation_nso_momentum_wheel_w_controller.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)
% $$$ 
% $$$ tempLibObj = 'Actuator Emulation/NSO momentum wheels';
% $$$ outputName = 'actuator_emulation_nso_momentum_wheels.eps';
% $$$ print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)

%%%%%%%%%% Complete model %%%%%%%%%%
tempLibObj = 'NSO';
outputName = 'complete_nso_model.eps';
print_model('',strcat(topLib,'/',tempLibObj),strcat(outputDir,'/',outputName),'portrait',1)


% Close top lib
close_system(topLib,0)
%%%%%%%%%%%%%%%%%%%%%% Print models %%%%%%%%%%%%%%%%%%%%
modelLib = '../models/';
model = 'basic_control_structure_kalman';
outputName = 'basis_model.eps';
print_model(modelLib,model,strcat(outputDir,'/',outputName),'portrait',1)

%Close remaining simulink system windows
bdclose('all');

fprintf(1,strcat('\n\t--= Simulink implementations printed and saved in :\n\t\t%s%s',outputDir,'\n'),current_dir,'\')
fprintf('\n')