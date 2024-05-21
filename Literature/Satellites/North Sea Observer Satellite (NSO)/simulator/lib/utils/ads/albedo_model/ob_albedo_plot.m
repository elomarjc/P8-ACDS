close all
clear all

%loading data file
ad = load('epr_data/2003/ga030101-031231.mat');

meas_deg = -89.5:1:89.5;
lat_mean = mean(ad.data');
lat_std = std(ad.data');

% polynomial fitting of data
p_order = 2;
[p_coef, p_structure] = polyfit(meas_deg,lat_mean,p_order);

%error
[p_val, p_err] = polyval(p_coef,meas_deg,p_structure);
rms_err = sqrt((p_err * p_err')/length(p_err));

%plotting data
fontsize = 10;
figure(1);
plot(meas_deg, lat_mean*100,'.',meas_deg, lat_std*100, ':', meas_deg, p_val*100, 'k');
title('Fitting of polynomial to latitude meaned albedo', 'FontSize', fontsize)
legend('Mean albedo','Standard deviation',[num2str(p_order) '. order polynomial fit'], 1);
ylabel('Reflectivity [%]', 'FontSize', fontsize);
xlabel('Latitude [\circ]', 'FontSize', fontsize);
text(-28,75,['rms_{err} on fit = ',num2str(rms_err*100), '%'])
axis([-90 90 0 100]);
set(gca,'XTick',[-90 -60 -30 0 30 60 90]);
grid;
shg