% Automated ADS test file
%
%
% Simulation data structure
% Description:                      Position        Unit
%--------------------------------------------------------
% Center of Mass                    (1:3)           [m]
% Mass of Satellite                 (4)             [kg]
% Inertia                           (5:7)           [kg*m]
% Q_S_C (q1,q2,q3,q4)               (8:11)          [-]
% Initial Attitude (q1,q2,q3,q4)    (12:15)         [-]
% Initial Angular Velocity          (16:18)         [rad/s]
% Initial JD (Real)                 (19)            [JD]
% Initial JD (SW)                   (20)            [JD]
% Sample Frequencies (data)         (21)            [Hz]
% Sample Frequencies (SW)           (22)            [Hz]
% Sample Frequencies (filters)      (23)            [Hz]
% Enable Permanent Magnet (1=true)  (24)            [-]
% SVD Observation Weights           (25:26)         [-]
% Magnetic Dipole Moment            (27)            [A*m^2]
%
%
% Sensor data structure
% Description:                  Default value       Unit
%--------------------------------------------------------
% Magnetometer:
% Sampling Frequency            1                   [Hz]
% Noise (standard deviation)    3                   [deg]
% Bias (x,y,z)                  (2000,4000,0)       [nT]
% Placement error (heading)     1                   [deg]
% placement error (tilt)        1                   [deg]

% Sunsensor:
% Sampling Frequency            1                   [Hz]
% Noise (standard deviation)    3.33                [deg]
% Bias (theta, phi)             (0,0)               [deg]
% Number of Samples (unused)    1                   [-]

% Gyroscope:
% Sampling Frequency            1                   [Hz]
% Noise (standard deviation)    0.2                 [deg]
% Bias (x,y,z)                  (-0.2,-0.2,-0.2)    [deg/s]
% Placement error (heading)     1                   [deg]
% placement error (tilt)        1                   [deg]

close all
clear
clc

%% Default Values

sensor_data=[1,3,[0,0,0],0,0,1,3.33,[0,0],1,1,0.2,[0,0,0],0,0];

% Simulation Time (start,end)
sim_time=[0,12000];

% Initial Julian Date for Simulation
init_JD_real=2455153.1954318+0.035;
init_JD_sw=2455153.1954318+0.035;

% Satellite parameters
CoM=[0.049073,0.048909,0.042976];
mass=0.957820;
inertia=[0.0017464,0.0022092,0.0022388];
q_s_c=[-0.020656,0.71468,-0.007319,0.6991];
init_att=[0 0 0 1];
init_w=[0.02,0.02,0.02];

% Sampling Frequencies (data,SW,filters)
f_sample=[1,1,1];

% Enable Permanent magnet (1=true,0=false)
enable_pmag=1;

% Magnetic Dipole Moment
mag_dipole=0.0030;

% SVD Weights for Observation Vectors (sun,mag)
svd_obs_weigth=[1/0.0034 1/0.0027];


%% Initilize
fprintf('Initializing... \n');
h_waitbar=waitbar(0,'Initializing...','name','Automated ADS simulation');

% Move to Test Folder
cd ../..;

% General
test_names={'UKF perfect conditions','UKF realistic bias mag',...
    'UKF realistic bias gyro','UKF realistic bias combo','UKF large bias mag',...
    'UKF large bias gyro','UKF large bias combo','UKF realistic inertia displacement',...
    'UKF unrealistic inertia displacement','UKF sensor displacement'};
save_names={'ukf_perf_cond','ukf_real_bias_mag','ukf_real_bias_gyro',...
    'ukf_real_bias_combo','ukf_large_bias_mag','ukf_large_bais_gyro',...
    'ukf_large_bias_combo','ukf_real_inertia','ukf_unreal_inertia','ukf_sensor_displacement'};

bar_steps = length(test_names)*4+3;
progress_count=1;
sim_access=0;
save_dir='test/ads/';
test_dir='test/ads/';
test_name='ads_test';


sat_par=[CoM,mass,inertia,q_s_c,init_att,init_w];
init_JD=[init_JD_real,init_JD_sw];
other_par=[f_sample,enable_pmag,svd_obs_weigth,mag_dipole];

sim_par=[sat_par,init_JD,other_par];

%% Check initial conditions
fprintf('Checking initial conditions... \n');
progress_count=progress_count+1;
waitbar(progress_count/bar_steps,h_waitbar,'Checking initial conditions...');

if length(test_names)==length(save_names)
    sim_check=1;
else
    sim_check=0;
end

%% Run simulation
if sim_check==1
    for i=1:length(test_names)
        fprintf('%s: Simulating... \n',char(test_names(i)));
        progress_count=progress_count+1;
        waitbar(progress_count/bar_steps,h_waitbar,[char(test_names(i)),': Simulating...']);
        
        switch i
            case 1 % Perfect
                sensor_data=[1,3,[0,0,0],0,0,1,3.33,[0,0],1,1,0.2,[0,0,0],0,0];
                sim_access=0;
            case 2 % Bias Mag (realistic)
                sensor_data=[1,3,[5000,1000,-3000],0,0,1,3.33,[0,0],1,1,0.2,[0,0,0],0,0];
                sim_access=1;
            case 3 % Bias Gyro (realistic)
                sensor_data=[1,3,[0,0,0],0,0,1,3.33,[0,0],1,1,0.2,[0.2,0.2,0.2],0,0];
                sim_access=0;
            case 4 % Bias Combo (realistic)
                sensor_data=[1,3,[5000,1000,-3000],0,0,1,3.33,[0,0],1,1,0.2,[0.2,0.2,0.2],0,0];
                sim_access=0;
            case 5 % Bias Mag (large)
                sensor_data=[1,3,[50000,10000,-30000],0,0,1,3.33,[0,0],1,1,0.2,[0,0,0],0,0];
                sim_access=0;
            case 6 % Bias Gyro (large)
                sensor_data=[1,3,[0,0,0],0,0,1,3.33,[0,0],1,1,0.2,[2.0,2.0,2.0],0,0];
                sim_access=0;
            case 7 % Bias Combo (large)
                sensor_data=[1,3,[50000,10000,-30000],0,0,1,3.33,[0,0],1,1,0.2,[2.0,2.0,2.0],0,0];
                sim_access=0;
            case 8 % Bias Combo (realistic) and Inertia (large)
                sensor_data=[1,3,[5000,1000,-3000],0,0,1,3.33,[0,0],1,1,0.2,[0.2,0.2,0.2],0,0];
                sim_par(1:3)=[0.049288,0.049161,0.07497]; % CoM
                sim_par(4)=1.25; %mass
                sim_par(5:7)=[0.0016986,0.002504,0.002545]; % inertia
                sim_par(8:11)=[-0.00033732,0.69934,-0.021604,0.71446]; % q_s_c
                sim_access=0;
            case 9 % Bias Combo (realistic) and Inertia (realistic)
                sensor_data=[1,3,[5000,1000,-3000],0,0,1,3.33,[0,0],1,1,0.2,[0.2,0.2,0.2],0,0];
                sim_par(1:3)=[0.0493,0.049176,0.034905]; % CoM
                sim_par(4)=1.27; %mass
                sim_par(5:7)=[0.0022634,0.0027338,0.0027634]; % inertia
                sim_par(8:11)=[-0.029144,0.71959,-0.003633,0.69378]; % q_s_c
                sim_access=0;
            case 10 % Bias Combo (realistic) and sensor displacement
                sensor_data=[1,3,[5000,1000,-3000],1,1,1,3.33,[0,0],1,1,0.2,[0.2,0.2,0.2],1,1];
                sim_par(1:3)=[0.049073,0.048909,0.042976]; % CoM
                sim_par(4)=0.957820; %mass
                sim_par(5:7)=[0.0017464,0.0022092,0.0022388]; % inertia
                sim_par(8:11)=[-0.020656,0.71468,-0.007319,0.6991]; % q_s_c
                sim_access=0;
            otherwise
                fprintf('%s is not a valid test! \n',char(test_names(i)));
                sim_access=0;
        end
        
        % Simulate
        if sim_access == 1
            sim([test_dir,test_name],sim_time);
        

            %% Store data
            fprintf('%s: Storing data... \n',char(test_names(i)));
            progress_count=progress_count+1;
            waitbar(progress_count/bar_steps,h_waitbar,[char(test_names(i)),': Storing data...']);
                
            data(:,63:65)=[tictoc_ukf(1,:)',tictoc_ukf_bias(1,:)',tictoc_svd(1,:)'];
            
            if exist([save_dir,char(save_names(i)),'.mat']) == 2
                for j=1:10
                    if exist([save_dir,char(save_names(i)),num2str(j),'.mat']) == 2
                        if j==10
                            save([save_dir,char(save_names(i)),num2str(j),'.mat'],'data','sensor_data','sim_par');
                        end
                    else
                        save([save_dir,char(save_names(i)),num2str(j),'.mat'],'data','sensor_data','sim_par');
                        break;
                    end
                end
            else
                save([save_dir,char(save_names(i)),'.mat'],'data','sensor_data','sim_par');
            end
            
            % Data structure
            %
            % Slot:                 Description:                    Unit:
            % Real
            %--------------------------------------------------------------
            % t(1)                  Time                            [s]
            % att(2:4)              Real attitude                   [deg]
            % w(5:7)                Real angular velocity           [deg/s]
            % UKF Simple
            %--------------------------------------------------------------
            % att_ukf(8:10)         Estimated attitude UKF          [deg]
            % w_ukf(11:13)          Estimated angular velocity      [deg/s]
            % r_sun(14:16)          Sun vector residual             [-]
            % r_mag(17:19)          Magnetic field vector residual  [-]
            % r_w(20:22)            Angular velocity residual       [deg/s]
            % eia(23:25)            Error in attitude               [deg]
            % eiw(26:28)            Error in angular velocity       [deg/s]
            % UKF Bias
            %--------------------------------------------------------------
            % att_ukf(29:31)        Estimated attitude UKF          [deg]
            % w_ukf(32:34)          Estimated angular velocity      [deg/s]
            % b_mag(35:37)          Bias magnetometer               [nT]
            % b_gyro(38:40)         Bias gyros                      [deg/s]
            % r_sun(41:43)          Sun vector residual             [-]
            % r_mag(44:46)          Magnetic field vector residual  [-]
            % r_w(47:49)            Angular velocity residual       [deg/s]
            % eia(50:52)            Error in attitude               [deg]
            % eiw(53:55)            Error in angular velocity       [deg/s]
            % SVD
            %--------------------------------------------------------------
            % att_svd(56:58)        Estimated attitude SVD          [deg]
            % e_att_svd(59:61)      Error in attitude for SVD       [deg]
            % svd_p_sum(62)         Sum of covariance matrix (SVD)  [-]
            % Other
            %--------------------------------------------------------------
            % tictoc_ukf(63)        Algorithm execution time (UKF)  [s]
            % tictoc_ukf_bias(64)   Algorithm execution time (Bias) [s]
            % tictoc_svd(65)
            % n_dist(66:68)         Environment disturbances        [Nm]
            % n_pmag(69:71)         Permanent magnet disturbance    [Nm]
            % n_con(72:74)          Control torque                  [Nm]
            % eclipse_real(75)      Eclipse (Real)                  [-]
            % eclipse_sw(76)        Eclipse (Software)              [-]



            %% Plot data
            fprintf('%s: Plotting data... \n',char(test_names(i)));
            progress_count=progress_count+1;
            waitbar(progress_count/bar_steps,h_waitbar,[char(test_names(i)),': Plotting data...'])

            pause(1);

            %% Save plots
            fprintf('%s: Saving plots... \n',char(test_names(i)));
            progress_count=progress_count+1;
            waitbar(progress_count/bar_steps,h_waitbar,[char(test_names(i)),': Saving plots...'])

            pause(1);
        else
            progress_count=progress_count+3;
        end

    end
end

%% Terminate
if sim_check == 1
    fprintf('Finished \n');
    progress_count=progress_count+1;
    waitbar(progress_count/bar_steps,h_waitbar,'Finished');
    cd(test_dir);
else
    fprintf('Simulation terminated... \n');
    waitbar(1,h_waitbar,'Simulation terminated...');
    cd(test_dir);
end

% Clean up
clear('init_JD_real','init_JD_sw','CoM','mass','inertia',...
    'q_s_c','init_att','init_w','f_sample','enable_pmag',...
    'svd_obs_weigth','other_par','init_JD','sat_par',...
    'mag_dipole','h_waitbar','test_names','save_names','bar_steps',...
    'progress_count','save_dir','test_dir','test_names','sim_check',...
    'sim_access','tout','test_name','i','ans','tictoc_ukf',...
    'tictoc_ukf_bias','tictoc_svd');
