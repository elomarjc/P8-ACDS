%% IGRF magnitude plot on IGRF map
% Clear command window
close all
clear
clc

% Run the 'igrf_test.mdl' file first with a total sample time of approximately 2 orbits
order = 13;
sim('igrf_test', 11650);

% Calculation of B field at 11 sample points and comparison with the program on the internet
% Program found on: http://www.ngdc.noaa.gov/geomagmodels/struts/calcIGRFWMM
% IGRF 10, 17/11-09 (spacecraft position from simulation has been used)
i = 7;
latitude_geodetic = sc_pos_phi(1+i*1000);
longitude = sc_pos_theta(1+i*1000);
Elevation = sc_pos_r(1+i*1000);

% Define differences between IGRF implementation and NOAA program
program(1) = 29381.2;
diff1 = num2str(B_field_tot(1+0*1000)-program(1),'%10.1f');
diff1 = sprintf(' %s', diff1);

program(2) = 42286.6;
diff2 = num2str(B_field_tot(1+1*1000)-program(2),'%10.1f');
diff2 = sprintf(' %s', diff2);

program(3) = 43875.5;
diff3 = num2str(B_field_tot(1+2*1000)-program(3),'%10.1f');
diff3 = sprintf(' %s', diff3);

program(4) = 22494.7;
diff4 = num2str(B_field_tot(1+3*1000)-program(4),'%10.1f');
diff4 = sprintf(' %s', diff4);
 
program(5) = 41675.2;
diff5 = num2str(B_field_tot(1+4*1000)-program(5),'%10.1f');
diff5 = sprintf(' %s', diff5);

program(6) = 35794.1;
diff6 = num2str(B_field_tot(1+5*1000)-program(6),'%10.1f');
diff6 = sprintf(' %s', diff6);

program(7) = 26766.5;
diff7 = num2str(B_field_tot(1+6*1000)-program(7),'%10.1f');
diff7 = sprintf(' %s', diff7);

program(8) = 41418.3;
diff8 = num2str(B_field_tot(1+7*1000)-program(8),'%10.1f');
diff8 = sprintf(' %s', diff8);

program(9) = 39325.2;
diff9 = num2str(B_field_tot(1+8*1000)-program(9),'%10.1f');
diff9 = sprintf(' %s', diff9);

program(10) = 25430.5;
diff10 = num2str(B_field_tot(1+9*1000)-program(10),'%10.1f');
diff10 = sprintf(' %s', diff10);

program(11) = 46095.4;
diff11 = num2str(B_field_tot(1+10*1000)-program(11),'%10.1f');
diff11 = sprintf(' %s', diff11);

program(12) = 25483.4;
diff12 = num2str(B_field_tot(1+11*1000)-program(12),'%10.1f');
diff12 = sprintf(' %s', diff12);

%% Plotting
figure
hold on;

% Size of font on x and y-label on plots
fontsize = 12;

% Create the 'background' axes (miller map projection)
ha = axesm('miller','PLineLocation',20); 

% Move the background axes to the bottom
uistack(ha,'bottom');      

% Load in a background image and display it using the correct colors 
I = imread('igrf_630km.jpg');              
hi = imagesc([-pi pi],[2.3 -2.3],I);
colormap gray

framem on; mlabel on; plabel on
setm(gca,'LabelFormat','signed','FontSize',fontsize)
axis normal
axis off

% Plot satellite trajectory and difference points
plotm(sc_pos_phi,sc_pos_theta,'.');

plotm(sc_pos_phi(1+0*1000),sc_pos_theta(1+0*1000),'Xr','MarkerSize',10,'LineWidth',2);
textm(sc_pos_phi(1+0*1000),sc_pos_theta(1+0*1000),diff1,'FontSize',fontsize,'FontWeight','demi','color','r')
plotm(sc_pos_phi(1+1*1000),sc_pos_theta(1+1*1000),'Xr','MarkerSize',10,'LineWidth',2);
textm(sc_pos_phi(1+1*1000),sc_pos_theta(1+1*1000),diff2,'FontSize',fontsize,'FontWeight','demi','color','r')
plotm(sc_pos_phi(1+2*1000),sc_pos_theta(1+2*1000),'Xr','MarkerSize',10,'LineWidth',2);
textm(sc_pos_phi(1+2*1000),sc_pos_theta(1+2*1000),diff3,'FontSize',fontsize,'FontWeight','demi','color','r')
plotm(sc_pos_phi(1+3*1000),sc_pos_theta(1+3*1000),'Xr','MarkerSize',10,'LineWidth',2);
textm(sc_pos_phi(1+3*1000),sc_pos_theta(1+3*1000),diff4,'FontSize',fontsize,'FontWeight','demi','color','r')
plotm(sc_pos_phi(1+4*1000),sc_pos_theta(1+4*1000),'Xr','MarkerSize',10,'LineWidth',2);
textm(sc_pos_phi(1+4*1000),sc_pos_theta(1+4*1000),diff5,'FontSize',fontsize,'FontWeight','demi','color','r')
plotm(sc_pos_phi(1+5*1000),sc_pos_theta(1+5*1000),'Xr','MarkerSize',10,'LineWidth',2);
textm(sc_pos_phi(1+5*1000),sc_pos_theta(1+5*1000),diff6,'FontSize',fontsize,'FontWeight','demi','color','r')
plotm(sc_pos_phi(1+6*1000),sc_pos_theta(1+6*1000),'Xr','MarkerSize',10,'LineWidth',2);
textm(sc_pos_phi(1+6*1000),sc_pos_theta(1+6*1000),diff7,'FontSize',fontsize,'FontWeight','demi','color','r')
plotm(sc_pos_phi(1+7*1000),sc_pos_theta(1+7*1000),'Xr','MarkerSize',10,'LineWidth',2);
textm(sc_pos_phi(1+7*1000),sc_pos_theta(1+7*1000),diff8,'FontSize',fontsize,'FontWeight','demi','color','r')
plotm(sc_pos_phi(1+8*1000),sc_pos_theta(1+8*1000),'Xr','MarkerSize',10,'LineWidth',2);
textm(sc_pos_phi(1+8*1000),sc_pos_theta(1+8*1000),diff9,'FontSize',fontsize,'FontWeight','demi','color','r')
plotm(sc_pos_phi(1+9*1000),sc_pos_theta(1+9*1000),'Xr','MarkerSize',10,'LineWidth',2);
textm(sc_pos_phi(1+9*1000),sc_pos_theta(1+9*1000),diff10,'FontSize',fontsize,'FontWeight','demi','color','r')
plotm(sc_pos_phi(1+10*1000),sc_pos_theta(1+10*1000),'Xr','MarkerSize',10,'LineWidth',2);
textm(sc_pos_phi(1+10*1000),sc_pos_theta(1+10*1000),diff11,'FontSize',fontsize,'FontWeight','demi','color','r')
plotm(sc_pos_phi(1+11*1000),sc_pos_theta(1+11*1000),'Xr','MarkerSize',10,'LineWidth',2);
textm(sc_pos_phi(1+11*1000),sc_pos_theta(1+11*1000),diff12,'FontSize',fontsize,'FontWeight','demi','color','r')

hl = legend('Total Magnetic Field Intensity Contours [nT]','Satellite trajectory','Difference (IGRF implementation - NOAA) [nT]','Location','Best');
set(hl,'FontSize',fontsize);

%% IGRF10 order comparison
simulation_time = 11650;

% Calculate magnetic fields for different IGRF orders
order = 13;
sim('igrf_test', simulation_time);
B_field_tot_13 = B_field_tot;

order = 10;
sim('igrf_test', simulation_time);
B_field_tot_10 = B_field_tot;

order = 9;
sim('igrf_test', simulation_time);
B_field_tot_9 = B_field_tot;

order = 8;
sim('igrf_test', simulation_time);
B_field_tot_8 = B_field_tot;

order = 7;
sim('igrf_test', simulation_time);
B_field_tot_7 = B_field_tot;

order = 6;
sim('igrf_test', simulation_time);
B_field_tot_6 = B_field_tot;

% Convert JD to orbits since start time
% JD 2455153.19542 is
% CE 2009 November 17 16:41:24.3 UT  Tuesday
orbits = (JD(:)-JD(1))*14.81823818;

figure
hold on
plot(orbits, B_field_tot_13-B_field_tot_10, 'y', 'LineWidth', 1.5)
plot(orbits, B_field_tot_13-B_field_tot_9, 'b', 'LineWidth', 1.5)
plot(orbits, B_field_tot_13-B_field_tot_8, 'g', 'LineWidth', 1.5)
plot(orbits, B_field_tot_13-B_field_tot_7, 'r', 'LineWidth', 1.5)
plot(orbits, B_field_tot_13-B_field_tot_6, 'm', 'LineWidth', 1.5)
grid
set(gca, 'FontSize', fontsize)
hl = legend('13. and 10. order', '13. and 9. order', '13. and 8. order', '13. and 7. order', '13. and 6. order', 'Location', 'Best');
set(hl, 'FontSize', 10);
xlabel('Orbit', 'FontSize', fontsize);
ylabel('Difference in total magnetic field intensity [nT]', 'FontSize', fontsize);
box on

% Standard deviation calculation
difference1 = B_field_tot_13 - B_field_tot_10;
mean1 = mean(difference1);
n = length(difference1);
for i = 1:n
    temp(i) = (difference1(i) - mean1)^2;
end
standard_deviation_13_10 = sqrt(sum(temp)/n);

difference2 = B_field_tot_13 - B_field_tot_9;
mean2 = mean(difference2);
n = length(difference2);
for i = 1:n
    temp(i) = (difference2(i) - mean2)^2;
end
standard_deviation_13_9 = sqrt(sum(temp)/n);

difference3 = B_field_tot_13 - B_field_tot_8;
mean3 = mean(difference3);
n = length(difference3);
for i = 1:n
    temp(i) = (difference3(i) - mean3)^2;
end
standard_deviation_13_8 = sqrt(sum(temp)/n);

difference4 = B_field_tot_13 - B_field_tot_7;
mean4 = mean(difference4);
n = length(difference4);
for i = 1:n
    temp(i) = (difference4(i) - mean4)^2;
end
standard_deviation_13_7 = sqrt(sum(temp)/n);

difference5 = B_field_tot_13 - B_field_tot_6;
mean5 = mean(difference5);
n = length(difference5);
for i = 1:n
    temp(i) = (difference5(i) - mean5)^2;
end
standard_deviation_13_6 = sqrt(sum(temp)/n);

%% Noise characteristics (difference between 8th and 13th order models)
simulation_time = 4 * 11650;

order = 13;
sim('igrf_test', simulation_time);
B_field_tot_13 = B_field_tot;

order = 8;
sim('igrf_test', simulation_time);
B_field_tot_8 = B_field_tot;

orbits = (JD(:)-JD(1)) * 14.81823818;

difference = B_field_tot_13 - B_field_tot_8;
figure
hold on
plot(orbits, difference, 'k', 'LineWidth', 1.5)
grid
set(gca, 'FontSize', fontsize)
xlabel('Orbit', 'FontSize', fontsize);
ylabel('Difference in total magnetic field intensity [nT]', 'FontSize', fontsize);
xlim([0 orbits(simulation_time)]);
ylim([-110 110]);
box on

% Standard deviation calculation
mean = mean(difference);
n = length(difference);
for i = 1:n
    temp(i) = (difference(i) - mean)^2;
end
standard_deviation = sqrt(sum(temp)/n);
max_difference = max(difference);
min_difference = min(difference);
