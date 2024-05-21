%This file is used for generating data for validation of the IGRF model
%under different conditions:
%
%   Condition 1:    The Ørsted position data is used with the IGRF model and compared
%                   to the Ørsted mag. data.
%   Condition 2:    The SPG4 model (with a appropriate TLE)is used with the IGRF model. 
%                   The results is compared to Ørsted mag. data 
%   Condition 3:    The SGP4 model is used with the IGRF model and compared to the
%                   Ørsted mag. data. The SGP4 model is disturbed by varying
%                   time-delays.
%
%   To get it working:
%   AAUSATII toolbox must be available (added to the path)
%   oersted_data.mat must be available (see below)

%% initialization
load oersted_data.mat % Ørsted pos. data and Ørsted mag. data.

fontsize=12;%size of font on x and y-label on plots
delay = 0;
modelorder=7;
simtime = length(oersted_data.time.sample_number);

%generate time data used for plotting
days = oersted_data.time.day(end) - oersted_data.time.day(1);%number of days in sim.data
start_sec = oersted_data.time.second(1); %start second of the first day of sim.data
end_sec = oersted_data.time.second(end); %end second of the last day  of sim.data
time_days = days + (end_sec - start_sec)/86400; %(86400=24*60*60)
total_sim_time_days = linspace(0, time_days, simtime+1);

%% prints the Simulink model (mag_field_model_test.mdl) to eps
open mag_field_model_test.mdl
print -smag_field_model_test -deps2 fig/igrf_val_simulink_rals.eps


%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation starts
sim('mag_field_model_test', simtime);
%%%%%%%%%%%%%%%%%%%%%%%%%


%Determine the maximal value for the for the magnitude plots in condition 1
%and 2 and convert to mGauss
mag_axis_max=1000*max(max(max(mag_error_orsted_mag_orsted_pos_igrf)),max(max(mag_error_sgp4_igrf_orsted_mag)));
mag_axis_min=1000*min(min(min(mag_error_orsted_mag_orsted_pos_igrf)),max(max(mag_error_sgp4_igrf_orsted_mag)));

%Determine the maximal value for the for the phase plots in condition 1
%and 2
phase_axis_max=max(max(phase_error_orsted_mag_orsted_pos_igrf),max(phase_error_sgp4_igrf_orsted_mag));


%% Condition 1
figure
hold on;
plot(total_sim_time_days,1000*mag_error_orsted_mag_orsted_pos_igrf);%plot in mGauss
ylabel('Error [mGauss]','fontsize',fontsize);
xlabel('Time [Day]','fontsize',fontsize);
axis([0 time_days mag_axis_min mag_axis_max]);
print -deps2c -r150 ./fig/mag_error_orsted_mag_orsted_pos_igrf_rals.eps

figure
hold on;
plot(total_sim_time_days,phase_error_orsted_mag_orsted_pos_igrf);
ylabel('Error [\circ]','fontsize',fontsize);
xlabel('Time [Day]','fontsize',fontsize);
axis([0 time_days 0 phase_axis_max]);
print -deps2c -r150 ./fig/phase_error_orsted_mag_orsted_pos_igrf_rals.eps


        Angular_error = phase_error_orsted_mag_orsted_pos_igrf;
        
        %%%%    MAGNITUDE ERROR:  %%%%
        ErrorfsMean = mean(mag_error_orsted_mag_orsted_pos_igrf);%RMS Value
        ErrorfsMx = max(mag_error_orsted_mag_orsted_pos_igrf);%Maximum value

        %%%%    ANGULAR ERROR:  %%%%
        ErroranMean = mean(Angular_error);
        ErroranMx = max(Angular_error);

        %%%%    RESULTS:  %%%%
        mag_error_orsted_mag_orsted_pos_igrf_MAX = 1000*ErrorfsMx; %Maximum value mGauss
        mag_error_orsted_mag_orsted_pos_igrf_RMS = 1000*ErrorfsMean; %RMS Value mGauss
        phase_error_orsted_mag_orsted_pos_igrf_MAX = ErroranMx; %Maximum value
        phase_error_orsted_mag_orsted_pos_igrf_RMS = ErroranMean; %RMS Value

%% Condition 2
figure
hold on;
% p_order=4;
% [coef,s]=polyfit(total_sim_time_days',mag_error_sgp4_igrf_orsted_mag*1000,p_order);
% plot(total_sim_time_days,mag_error_sgp4_igrf_orsted_mag*1000,'b',...
%     total_sim_time_days,polyval(coef,total_sim_time_days),'r');% in km
plot(total_sim_time_days,1000*mag_error_sgp4_igrf_orsted_mag);%plot in mGauss
ylabel('Error [mGauss]','fontsize',fontsize);
xlabel('Time [Day]','fontsize',fontsize);
axis([0 time_days mag_axis_min mag_axis_max]);
%legend('Error, normed magnetic field vectors','4th order polyfit',2)
print -deps2c -r150 ./fig/mag_error_sgp4_igrf_orsted_mag_rals.eps

figure
hold on;
% p_order=4;
% [coef,s]=polyfit(total_sim_time_days',phase_error_sgp4_igrf_orsted_mag,p_order);
% plot(total_sim_time_days,phase_error_sgp4_igrf_orsted_mag,'b',...
%     total_sim_time_days,polyval(coef,total_sim_time_days),'r');% in km
plot(total_sim_time_days,phase_error_sgp4_igrf_orsted_mag);
ylabel('Error [\circ]','fontsize',fontsize);
xlabel('Time [Day]','fontsize',fontsize);
axis([0 time_days 0 phase_axis_max]);
%legend('Angular error','4th order polyfit',2)
print -deps2c -r150 ./fig/phase_error_sgp4_igrf_orsted_mag_rals.eps

        
        %%%%    MAGNITUDE ERROR:  %%%%
        ErrorfsMean = mean(mag_error_sgp4_igrf_orsted_mag);%RMS Value
        ErrorfsMx = max(mag_error_sgp4_igrf_orsted_mag);%Maximum value

        %%%%    ANGULAR ERROR:  %%%%
        ErroranMean = mean(phase_error_sgp4_igrf_orsted_mag);
        ErroranMx = max(phase_error_sgp4_igrf_orsted_mag);

        %%%%    RESULTS:  %%%%
        mag_error_sgp4_igrf_orsted_mag_MAX = 1000*ErrorfsMx; %Maximum value mGauss
        mag_error_sgp4_igrf_orsted_mag_RMS = 1000*ErrorfsMean; %RMS Value mGauss
        phase_error_sgp4_igrf_orsted_mag_MAX = ErroranMx; %Maximum value
        phase_error_sgp4_igrf_orsted_mag_RMS = ErroranMean; %RMS Value

% Calculate RMS and MAX error pr. day
samples_pr_day = round(sample_size/time_days);
% 
% for i = 1:1:7
%     
%         %%%%    MAGNITUDE ERROR:  %%%%
%         tempfs = field_strength_error.^2; 
%         tempfs2 = tempfs(:,1)+tempfs(:,2)+tempfs(:,3);
%         ErrorfsMean = mean(sqrt(tempfs2(1+samples_pr_day*(i-1):samples_pr_day*i)));%RMS Value
%         ErrorfsMx = max(max(field_strength_error(1+samples_pr_day*(i-1):samples_pr_day*i,:)));%Maximum value
% 
%         %%%%    ANGULAR ERROR:  %%%%
%         ErroranMean = sqrt(mean(Angular_error(1+samples_pr_day*(i-1):samples_pr_day*i).^2));
%         ErroranMx = max(Angular_error(1+samples_pr_day*(i-1):samples_pr_day*i,:));
% 
%         %%%%    RESULTS:  %%%%
%         ErrorArray_days(i,1) = i;
%         ErrorArray_days(i,2) = 1000*ErrorfsMx;     %Magnitude Maximum value mGauss
%         ErrorArray_days(i,3) = 1000*ErrorfsMean;   %Magnitude RMS Value mGauss
%         ErrorArray_days(i,4) = ErroranMx;     %Phase Maximum value
%         ErrorArray_days(i,5) = ErroranMean;   %Phase RMS Value
% end
% ErrorArray_days(:,3)/1000
% ErrorArray_days(:,2)/1000
% ErrorArray_days(:,5)
% ErrorArray_days(:,4)


%% Condition 3
sim_delay = 0;
if sim_delay == 1
    for i=0.5:0.5:5
        delay = i;

simtime = length(oersted_data.time.sample_number);

        %%%%%%%%%%%%%%%%%%%%%%%%%
        %Simulation starts
        sim('mag_field_model_test', simtime);
        %%%%%%%%%%%%%%%%%%%%%%%%%
        field_strength_error = mag_error_sgp4_igrf_orsted_mag;
        Angular_error = phase_error_sgp4_igrf_orsted_mag;

        %%%%    MAGNITUDE ERROR:  %%%%

        ErrorfsMean = mean(mag_error_sgp4_igrf_orsted_mag);%RMS Value
        ErrorfsMx = max(mag_error_sgp4_igrf_orsted_mag);%Maximum value

        %%%%    ANGULAR ERROR:  %%%%

        ErroranMean = mean(phase_error_sgp4_igrf_orsted_mag);
        ErroranMx = max(phase_error_sgp4_igrf_orsted_mag);


        %%%%    RESULTS:  %%%%
        Index = 2*i; %only to be used with i=0.5:0.5:5 ;-)

        ErrorArray_delay(Index,1) = i;
        ErrorArray_delay(Index,2) = ErrorfsMx; %Maximum value
        ErrorArray_delay(Index,3) = ErrorfsMean; %RMS Value 
        ErrorArray_delay(Index,4) = ErroranMx; %Maximum value
        ErrorArray_delay(Index,5) = ErroranMean %RMS Value

        %plot(field_strength_error)
        %axis([0 sample_size min(min(field_strength_error)) max(max(field_strength_error))]);
        %xlabel('Samples');
        %ylabel('Error [mGauss]');
        %Title = ['Order' int2str(order)]

    end

end
        %save ErrorArray_sgp4_delay_igrf.mat ErrorArray_delay %The simulation time is approx. ½ hour on this pc for the entire script

        
asdf = ErrorArray_delay'*1000
for idx=1:size(asdf,1), disp(['MEasd ' sprintf(' & %.2f ',asdf(idx,:)) '\\']),end

mag_error_orsted_mag_orsted_pos_igrf_MAX %Maximum value mGauss
mag_error_orsted_mag_orsted_pos_igrf_RMS %RMS Value mGauss
phase_error_orsted_mag_orsted_pos_igrf_MAX %Maximum value
phase_error_orsted_mag_orsted_pos_igrf_RMS %RMS Value


mag_error_sgp4_igrf_orsted_mag_MAX %Maximum value mGauss
mag_error_sgp4_igrf_orsted_mag_RMS %RMS Value mGauss
phase_error_sgp4_igrf_orsted_mag_MAX %Maximum value
phase_error_sgp4_igrf_orsted_mag_RMS %RMS Value


%% upload eps-files
disp('uploading files...')
!pscp.exe -batch fig/*.* AAU:cvs/adcs_common/appendix/fig/

