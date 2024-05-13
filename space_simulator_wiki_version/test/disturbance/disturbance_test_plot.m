%% Disturbance test plot
close all
clear
clc

fontsize=12;

%simulate one orbit
sim('disturbance_test', 5830*10);

orbits=(JD(:)-JD(1))*14.81823818;
time_sec=(JD(:)-JD(1))*86400;

%% Figure 1
figure

%N_earth 
h = subplot(4,1,1);
set(h,'position',[0.1 0.05+0.03+0.20+0.20+0.03+0.20+0.03 0.85 0.20]); 
plot(orbits,sqrt(N_earth(:,1).^2+N_earth(:,2).^2+N_earth(:,3).^2))
set(gca,'xtick',[])
ylabel('Torque [Nm]','FontSize',fontsize);
legend('Source: Earth');
max_gravity_gradient=max(sqrt(N_earth(:,1).^2+N_earth(:,2).^2+N_earth(:,3).^2))
%title('Magnitude of Gravity Gradient Disturbance Torques','FontSize',fontsize);

%N_moon 
h=subplot(4,1,2);
set(h,'position',[0.1 0.05+0.03+0.20+0.20+0.03 0.85 0.20]); 
plot(orbits,sqrt(N_moon(:,1).^2+N_moon(:,2).^2+N_moon(:,3).^2))
set(gca,'xtick',[])
ylabel('Torque [Nm]','FontSize',fontsize);
legend('Source: Moon');
max_moon=max(sqrt(N_moon(:,1).^2+N_moon(:,2).^2+N_moon(:,3).^2))
%title('Magnitmax_moon=max(sqrt(N_moon(:,1).^2+N_moon(:,2).^2+N_moon(:,3).^2))ude of Moon Disturbance Torques','FontSize',fontsize);

%N_sun
h=subplot(4,1,3);
set(h,'position',[0.1 0.05+0.03+0.20 0.85 0.20]); 
plot(orbits,sqrt(N_sun(:,1).^2+N_sun(:,2).^2+N_sun(:,3).^2))
set(gca,'xtick',[])
ylabel('Torque [Nm]','FontSize',fontsize);
legend('Source: Sun');
max_sun=max(sqrt(N_sun(:,1).^2+N_sun(:,2).^2+N_sun(:,3).^2))
%title('Magnitude of Sun Disturbance Torques','FontSize',fontsize);

%N_zonal
h=subplot(4,1,4);
set(h,'position',[0.1 0.05 0.85 0.20]); 
plot(orbits,sqrt(N_zonal(:,1).^2+N_zonal(:,2).^2+N_zonal(:,3).^2))
xlabel('Orbit [-]','FontSize',fontsize);
ylabel('Torque [Nm]','FontSize',fontsize);
legend('Source: Zonal harmonics');
max_zonal=max(sqrt(N_zonal(:,1).^2+N_zonal(:,2).^2+N_zonal(:,3).^2))
%title('Magnitude of Zonal Harmonics Disturbance Torques','FontSize',fontsize);

%% Figure 2
figure

%N_mag_res
h=subplot(3,1,1);
set(h,'position',[0.1 0.065+0.25+0.05+0.25+0.05 0.85 0.26]);
plot(orbits,sqrt(N_mag_res(:,1).^2+N_mag_res(:,2).^2+N_mag_res(:,3).^2))
set(gca,'xtick',[])
ylabel('Torque [Nm]','FontSize',fontsize);
legend('Source: Magnetic residual');
max_mag_res=max(sqrt(N_mag_res(:,1).^2+N_mag_res(:,2).^2+N_mag_res(:,3).^2))
%title('Magnitude of Magnetic Residual Disturbance Torques','FontSize',fontsize);

%N_atmospheric
h=subplot(3,1,2);
set(h,'position',[0.1 0.065+0.25+0.05 0.85 0.26]);
plot(orbits,sqrt(N_atmospheric(:,1).^2+N_atmospheric(:,2).^2+N_atmospheric(:,3).^2))
set(gca,'xtick',[])
ylabel('Torque [Nm]','FontSize',fontsize);
legend('Source: Atmospheric drag');
max_drag=max(sqrt(N_atmospheric(:,1).^2+N_atmospheric(:,2).^2+N_atmospheric(:,3).^2))
%title('Magnitude of Atmospheric Drag Disturbance Torques','FontSize',fontsize);

%N_radiation
h=subplot(3,1,3);
set(h,'position',[0.1 0.065 0.85 0.26]);
hold on
plot(orbits,0.5*10^(-9)*eclipse,'r')
plot(orbits,sqrt(N_radiation(:,1).^2+N_radiation(:,2).^2+N_radiation(:,3).^2))
xlabel('Orbit [-]','FontSize',fontsize);
ylim([0 1.4*10^(-9)])
ylabel('Torque [Nm]','FontSize',fontsize);
legend('Eclipse indication (high=eclipse)','Source: Solar radiation');
box on
max_radiation=max(sqrt(N_radiation(:,1).^2+N_radiation(:,2).^2+N_radiation(:,3).^2))
%title('Magnitude of Solar Radiation Disturbance Torques','FontSize',fontsize);

%% Figure 3
figure

%N_totol
plot(orbits,disturbance_tot)
xlabel('Orbit [-]','FontSize',fontsize);
ylabel('Torque [Nm]','FontSize',fontsize);
%legend('Source: Total disturbance torque');
max_total=max(disturbance_tot)
%title('Magnitude of All Disturbance Torques','FontSize',fontsize);

%% Figure 4
figure

%Angular velocity
h=subplot(3,1,1);
set(h,'position',[0.1 0.065+0.25+0.065+0.25+0.065 0.85 0.26]);
plot(orbits,180/pi*omega_s(:,1),orbits,180/pi*omega_s(:,2),orbits,180/pi*omega_s(:,3))
xlabel('Orbit [-]','FontSize',fontsize);
ylabel('Angular velocity [deg/s]','FontSize',fontsize);
legend('X','Y','Z');
%title('Angular Velocity in SBRF','FontSize',fontsize);

%Attitude
h=subplot(3,1,2);
set(h,'position',[0.1 0.065+0.25+0.065 0.85 0.26]);
plot(orbits,180/pi*attitude(:,1),'X',orbits,180/pi*attitude(:,2),'X',orbits,180/pi*attitude(:,3),'X','MarkerSize',2)
xlabel('Orbit [-]','FontSize',fontsize);
ylim([-180 180]);
ylabel('Attitude [deg]','FontSize',fontsize);
legend('X','Y','Z');
%title('Attitude ECI to SBRF','FontSize',fontsize);

%Satellite trajectory
h=subplot(3,1,3);
set(h,'position',[0.1 0.065 0.85 0.26]);
hold on
rgb = imread('earth_map.jpg'); %picture generated with orbitron. The white star correponds to the TLE timestamp
image([-180 180],[90 -90],rgb);
set(gca,'YDir','normal')
set(gca,'ytick',[-90 -60 -30 0 30 60 90])
set(gca,'xtick',[-180 -150 -120 -90 -60 -30 0 30 60 90 120 150 180])
xlim([-180 180]);
ylim([-90 90]);
plot(sc_longitude,sc_latitude,'.')
xlabel('Longitude [deg]','FontSize',fontsize);
ylabel('Latitude [deg]','FontSize',fontsize);
legend('Satellite trajectory');
%title('Satellite Trajectory','FontSize',fontsize);


