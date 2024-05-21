%   To get it working:s
%   AAUSATII toolbox must be available (added to the path)
%   oersted_data.mat must be available (see below)

disp('Remember to create a library named "fig"');

%% print simulink model to eps
open_system('kepler_sgp4_test')
print -s'kepler_sgp4_test' -deps2 fig/orbit_precision_simulink_rals.eps
close_system('kepler_sgp4_test')

%% init
load oersted_data.mat % Ørsted pos. data. 


%generate time data used for plotting
days = oersted_data.time.day(end) - oersted_data.time.day(1);%number of days in sim.data
start_sec = oersted_data.time.second(1); %start second of the first day of sim.data
end_sec = oersted_data.time.second(end); %end second of the last day  of sim.data
time_days = days + (end_sec - start_sec)/86400; %(86400=24*60*60)
total_sim_time_days = linspace(0, time_days, length(oersted_data.time.sample_number)+1);
sample_size=length(oersted_data.time.sample_number);

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation starts
sim('kepler_sgp4_test', length(oersted_data.time.sample_number));
%%%%%%%%%%%%%%%%%%%%%%%%%

fontsize=12;%size of font on x and y-label on plots

%Determine the maximal value for the for the magnitude plots in condition 1
%and 2
mag_axis_max=max(max(max(Error_SGP4_magnitude_workspace)),max(max(Error_kepler_magnitude_workspace)))/1000;%in km

%Determine the maximal value for the for the phase plots in condition 1
%and 2
phase_axis_max=max(max(Error_kepler_angle_workspace),max(Error_kepler_angle_workspace));

%% plot Error_SGP4_magnitude
figure
hold on;
p_order=4;
[coef,s]=polyfit(total_sim_time_days',Error_SGP4_magnitude_workspace/1000,p_order);
plot(total_sim_time_days,Error_SGP4_magnitude_workspace/1000,'b',...
    total_sim_time_days,polyval(coef,total_sim_time_days),'r');% in km
ylabel('Error [km]','fontsize',fontsize);
xlabel('Time [Day]','fontsize',fontsize);
axis([0 time_days 0 mag_axis_max]);
legend('Norm of position error','4th order polyfit',2)
print -depsc2 -r150 './fig/error_sgp4_mag_rals.eps'

%% plot Error_SGP4_angle
figure
hold on;
p_order=4;
[coef,s]=polyfit(total_sim_time_days',Error_SGP4_angle_workspace,p_order);
plot(total_sim_time_days,Error_SGP4_angle_workspace,'b',...
    total_sim_time_days,polyval(coef,total_sim_time_days),'r');% in km
ylabel('Error [\circ]','fontsize',fontsize);
xlabel('Time [Day]','fontsize',fontsize);
axis([0 time_days 0 phase_axis_max]);
legend('Angular error','4th order polyfit',2)
print -depsc2 -r150 './fig/error_sgp4_angle_rals.eps'

%% Error_kepler_magnitude
figure
hold on;
p_order=4;
[coef,s]=polyfit(total_sim_time_days',Error_kepler_magnitude_workspace/1000,p_order);
plot(total_sim_time_days,Error_kepler_magnitude_workspace/1000,'b',...
    total_sim_time_days,polyval(coef,total_sim_time_days),'r');% in km
ylabel('Error [km]','fontsize',fontsize);
xlabel('Time [Day]','fontsize',fontsize);
axis([0 time_days 0 mag_axis_max]);
legend('Norm of position error','4th order polyfit',2)
print -depsc2 -r150 './fig/error_kepler_mag_rals.eps'

%% Error_kepler_angle
figure
hold on;
p_order=4;
[coef,s]=polyfit(total_sim_time_days',Error_kepler_angle_workspace,p_order);
plot(total_sim_time_days,Error_kepler_angle_workspace,'b',...
    total_sim_time_days,polyval(coef,total_sim_time_days),'r');% in km
ylabel('Error [\circ]','fontsize',fontsize);
xlabel('Time [Day]','fontsize',fontsize);
axis([0 time_days 0 phase_axis_max]);
legend('Angular error','4th order polyfit',2)
print -depsc2 -r150 './fig/error_kepler_angle_rals.eps'


