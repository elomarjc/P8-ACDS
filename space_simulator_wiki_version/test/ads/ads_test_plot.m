    function ads_test_plot(test_name,check_data,varargin)
% Plot ADS test data
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

%% Initial
close all
clc

%% load data
load(test_name)
t_filters=[mean(data(:,63)),mean(data(:,64)),mean(data(:,65))];
load('ukf_real_bias_combo_no_qsc')
t_filters=[t_filters mean(data(:,64))];
fprintf('\nAlgorithm Execution Time - UKF: %4.2f [ms], UKF Bias: %4.2f [ms], SVD: %4.2f [ms], UKF Bias (no qsc): %4.2f [ms]\n',t_filters*1000);

load(test_name)

%% Plot True Attitude

plot_nr=0;

fig_size=[6,21]; % cm
figure('PaperUnits', 'centimeters','PaperPosition', [0.0,(29.2-fig_size(1)) fig_size(2) fig_size(1)]); %210 � 297 A4
set(gcf,'DefaultAxesFontSize',8,'DefaultAxesFontName', 'Times New Roman');
set(gcf,'DefaultAxesColorOrder',[0,0,0;0.6,0.6,0.6],'DefaultAxesLineStyleOrder','-|--');
hold on

temp=sign(data(1:12000,2:4));
par={'-','--','-'};
parc={[0,0,0],[0,0,0],[0.4,0.4,0.4]};
for k=1:3
    temp2=[];
    temp2=[temp2 0];
    for i=1:length(temp(:,k))-1
        if temp(i,k)>0 && temp(i+1,k)<0 && data(i,k+1)>150 && data(i+1,k+1)<-150
            temp2=[temp2 i];
        end
        if temp(i+1,k)>0 && temp(i,k)<0 && data(i+1,k+1)>150 && data(i,k+1)<-150
            temp2=[temp2 i];
        end
    end
    temp2=[temp2 12000];
    for j=1:length(temp2)-1
        plot(data(temp2(j)+1:temp2(j+1),1),data(temp2(j)+1:temp2(j+1),k+1),'LineWidth',0.9,'LineStyle',par{k},'Color',parc{k});
    end
    legend('x-axis','Location','NorthEast');
end
set(gca,'YTick',(-180:60:180))
axis([0 6000 -180 180]);
xlabel('Time [s]')
ylabel('Attitude [deg]')
grid;
hold off

if nargin > 3
    if strcmp(varargin(1),{'plotname'}) == 1
        print('-dpdf',[varargin{2},num2str(plot_nr)])
    end
end

%fprintf('Legnth of temp2 is %f \n',length(temp2))
%fprintf('Legnth of temp2 is %s \n',num2str(temp2(1:10)))


%% Plot SVD and Sum of Covariance

load('ukf_real_bias_combo3')

s_b=1;
for ii=1:length(data(:,75))
    if data(ii,75)>0 && s_b==1
        eclip(1)=ii;
        s_b=0;
    end
    if data(ii,75)<1 && s_b==0
        eclip(2)=ii-1;
        break;
    end
end
eclip(1)
eclip(2)
stage(1)=max(max(abs(data(1:1000,59:61))));
stage(2)=max(max(abs(data(1000+1:2100,59:61))));
stage(3)=max(max(abs(data(2101:4095,59:61))));
stage(4)=max(max(abs(data(4096:6000,59:61))));

fprintf('\nMax deviations:\n');
fprintf('\nPlot:      Stage 1:         Stage 2:         Stage 3:         Stage 4:\n');
fprintf('\n  1%12.2f [deg]%12.2f [deg]%12.2f [deg]%12.2f [deg]\n',stage(1:4));

plot_nr=1;
fig_size=[13,21]; % cm
figure('PaperUnits', 'centimeters','PaperPosition', [0.0,(29.2-fig_size(1)) fig_size(2) fig_size(1)]); %210 � 297 A4
set(gcf,'DefaultAxesFontSize',8,'DefaultAxesFontName', 'Times New Roman');
set(gcf,'DefaultAxesColorOrder',[0,0,0;0.6,0.6,0.6],'DefaultAxesLineStyleOrder','-|--');
%set(gcf,'PaperSize',fig_size*2.54);
subplot(2,1,1);
plot(data(:,1),data(:,59:61),'LineWidth',0.9);
set(gca,'XtickLabel',[])
set(gca,'ytick',(-180:30:180))
ylabel('Error in Estimated Attitude [deg]')
legend('x','y','z')
axis([0 6000 -180 180]);
grid;

set(gcf,'DefaultAxesFontSize',8,'DefaultAxesFontName', 'Times New Roman');
subplot(2,1,2)
h=get(gca,'Position');
set(gca,'Position',[h(1) h(2)+0.1 h(3) h(4)])
plot(data(:,1),data(:,62),'LineWidth',0.9);
xlabel('Time [s]')
ylabel('Sum of the Covariance Matrix [-]')
axis([0 6000 0 0.5]);
grid;

if nargin > 3
    if strcmp(varargin(1),{'plotname'}) == 1
        print('-dpdf',[varargin{2},num2str(plot_nr)])
    end
end

%% Plot UKF and UKF Bias without Bias and with Realistic Bias (combo)

load(test_name)

s_b=1;
for ii=1:length(data(:,75))
    if data(ii,75)>0 && s_b==1
        eclip(1)=ii;
        s_b=0;
    end
    if data(ii,75)<1 && s_b==0
        eclip(2)=ii-1;
        break;
    end
end

stage(1)=max(max(abs(data(1:1000,23:25))));
stage(2)=max(max(abs(data(1000+1:2100,23:25))));
stage(3)=max(max(abs(data(2101:4095,23:25))));
stage(4)=max(max(abs(data(4096:6000,23:25))));
fprintf('\n  2%12.2f [deg]%12.2f [deg]%12.2f [deg]%12.2f [deg]\n',stage(1:4));

plot_nr=2;
fig_size=[28,21]; % cm
figure('PaperUnits', 'centimeters','PaperPosition', [0.0,(29.2-fig_size(1)) fig_size(2) fig_size(1)]); %210 � 297 A4
set(gcf,'DefaultAxesFontSize',8,'DefaultAxesFontName', 'Times New Roman');
set(gcf,'DefaultAxesColorOrder',[0,0,0;0.6,0.6,0.6],'DefaultAxesLineStyleOrder','-|--');
subplot(4,1,1);
plot(data(:,1),data(:,23:25),[eclip(1);eclip(2)],[-9 -9],'k-',[eclip(1);eclip(1)],[-9.5 -8.5],'k-',[eclip(2);eclip(2)],[-9.5 -8.5],'k-','LineWidth',0.9);
text(eclip(1)+100,-8,'Eclipse','FontSize',8,'FontName','Times New Roman');
set(gca,'XtickLabel',[])
set(gca,'ytick',(-10:5:10))
ylabel('Error in Estimated Attitude [deg]')
legend('x','y','z')
axis([0 6000 -10 10]);
grid;


stage(1)=max(max(abs(data(1:1000,50:52))));
stage(2)=max(max(abs(data(1000+1:2100,50:52))));
stage(3)=max(max(abs(data(2101:4095,50:52))));
stage(4)=max(max(abs(data(4096:6000,50:52))));
fprintf('\n  3%12.2f [deg]%12.2f [deg]%12.2f [deg]%12.2f [deg]\n',stage(1:4));

set(gcf,'DefaultAxesFontSize',8,'DefaultAxesFontName', 'Times New Roman');
subplot(4,1,2);
h=get(gca,'Position');
set(gca,'Position',[h(1) h(2)+0.04 h(3) h(4)])
plot(data(:,1),data(:,50:52),[eclip(1);eclip(2)],[-9 -9],'k-',[eclip(1);eclip(1)],[-9.5 -8.5],'k-',[eclip(2);eclip(2)],[-9.5 -8.5],'k-','LineWidth',0.9);
text(eclip(1)+100,-8,'Eclipse','FontSize',8,'FontName','Times New Roman');
set(gca,'XtickLabel',[])
set(gca,'ytick',(-10:5:10))
ylabel('Error in Estimated Attitude [deg]')
legend('x','y','z')
axis([0 6000 -10 10]);
grid;

load('ukf_real_bias_combo')

s_b=1;
for ii=1:length(data(:,75))
    if data(ii,75)>0 && s_b==1
        eclip(1)=ii;
        s_b=0;
    end
    if data(ii,75)<1 && s_b==0
        eclip(2)=ii-1;
        break;
    end
end

stage(1)=max(max(abs(data(1:1000,23:25))));
stage(2)=max(max(abs(data(1000+1:2100,23:25))));
stage(3)=max(max(abs(data(2101:4095,23:25))));
stage(4)=max(max(abs(data(4096:6000,23:25))));
fprintf('\n  4%12.2f [deg]%12.2f [deg]%12.2f [deg]%12.2f [deg]\n',stage(1:4));

set(gcf,'DefaultAxesFontSize',8,'DefaultAxesFontName', 'Times New Roman');
subplot(4,1,3);
h=get(gca,'Position');
set(gca,'Position',[h(1) h(2)+0.08 h(3) h(4)])
plot(data(:,1),data(:,23:25),[eclip(1);eclip(2)],[-108 -108],'k-',[eclip(1);eclip(1)],[-114 -102],'k-',[eclip(2);eclip(2)],[-114 -102],'k-','LineWidth',0.9);
text(eclip(1)+100,-96,'Eclipse','FontSize',8,'FontName','Times New Roman');
set(gca,'XtickLabel',[])
set(gca,'ytick',(-120:60:120))
ylabel('Error in Estimated Attitude [deg]')
legend('x','y','z')
axis([0 6000 -120 120]);
grid;

stage(1)=max(max(abs(data(1:1000,50:52))));
stage(2)=max(max(abs(data(1000+1:2100,50:52))));
stage(3)=max(max(abs(data(2101:4095,50:52))));
stage(4)=max(max(abs(data(4096:6000,50:52))));
fprintf('\n  5%12.2f [deg]%12.2f [deg]%12.2f [deg]%12.2f [deg]\n',stage(1:4));

set(gcf,'DefaultAxesFontSize',8,'DefaultAxesFontName', 'Times New Roman');
subplot(4,1,4);
h=get(gca,'Position');
set(gca,'Position',[h(1) h(2)+0.12 h(3) h(4)])
plot(data(:,1),data(:,50:52),[eclip(1);eclip(2)],[-9 -9],'k-',[eclip(1);eclip(1)],[-9.5 -8.5],'k-',[eclip(2);eclip(2)],[-9.5 -8.5],'k-','LineWidth',0.9);
text(eclip(1)+100,-8,'Eclipse','FontSize',8,'FontName','Times New Roman');
xlabel('Time [s]')
ylabel('Error in Estimated Attitude [deg]')
legend('x','y','z')
set(gca,'ytick',(-10:5:10))
axis([0 6000 -10 10]);
grid;

if nargin > 3
    if strcmp(varargin(1),{'plotname'}) == 1
        print('-dpdf',[varargin{2},num2str(plot_nr)])
    end
end

%% Plot UKF Bias with Realistic Inertia and UKF Bias with no qsc

load('ukf_real_inertia')

s_b=1;
for ii=1:length(data(:,75))
    if data(ii,75)>0 && s_b==1
        eclip(1)=ii;
        s_b=0;
    end
    if data(ii,75)<1 && s_b==0
        eclip(2)=ii-1;
        break;
    end
end

stage(1)=max(max(abs(data(1:1000,50:52))));
stage(2)=max(max(abs(data(1000+1:2100,50:52))));
stage(3)=max(max(abs(data(2101:4095,50:52))));
stage(4)=max(max(abs(data(4096:6000,50:52))));
fprintf('\n  6%12.2f [deg]%12.2f [deg]%12.2f [deg]%12.2f [deg]\n',stage(1:4));

plot_nr=3;
fig_size=[13,21]; % cm
figure('PaperUnits', 'centimeters','PaperPosition', [0.0,(29.2-fig_size(1)) fig_size(2) fig_size(1)]); %210 � 297 A4
set(gcf,'DefaultAxesColorOrder',[0,0,0;0.6,0.6,0.6],'DefaultAxesLineStyleOrder','-|--');
set(gcf,'DefaultAxesFontSize',8,'DefaultAxesFontName', 'Times New Roman');
subplot(2,1,1);
plot(data(:,1),data(:,50:52),[eclip(1);eclip(2)],[-9 -9],'k-',[eclip(1);eclip(1)],[-9.5 -8.5],'k-',[eclip(2);eclip(2)],[-9.5 -8.5],'k-','LineWidth',0.9);
text(eclip(1)+100,-8,'Eclipse','FontSize',8,'FontName','Times New Roman');
set(gca,'XtickLabel',[])
ylabel('Error in Estimated Attitude [deg]');
legend('x','y','z')
axis([0 6000 -10 10]);
grid;

load('ukf_real_bias_combo_no_qsc')

s_b=1;
for ii=1:length(data(:,75))
    if data(ii,75)>0 && s_b==1
        eclip(1)=ii;
        s_b=0;
    end
    if data(ii,75)<1 && s_b==0
        eclip(2)=ii-1;
        break;
    end
end

stage(1)=max(max(abs(data(1:1000,50:52))));
stage(2)=max(max(abs(data(1000+1:2100,50:52))));
stage(3)=max(max(abs(data(2101:4095,50:52))));
stage(4)=max(max(abs(data(4096:6000,50:52))));
fprintf('\n  7%12.2f [deg]%12.2f [deg]%12.2f [deg]%12.2f [deg]\n',stage(1:4));

set(gcf,'DefaultAxesFontSize',8,'DefaultAxesFontName', 'Times New Roman');
subplot(2,1,2);
h=get(gca,'Position');
set(gca,'Position',[h(1) h(2)+0.1 h(3) h(4)])
plot(data(:,1),data(:,50:52),[eclip(1);eclip(2)],[-9 -9],'k-',[eclip(1);eclip(1)],[-9.5 -8.5],'k-',[eclip(2);eclip(2)],[-9.5 -8.5],'k-','LineWidth',0.9);
text(eclip(1)+100,-8,'Eclipse','FontSize',8,'FontName','Times New Roman');
xlabel('Time [s]');
ylabel('Error in Estimated Attitude [deg]');
legend('x','y','z')
axis([0 6000 -10 10]);
grid;

if nargin > 3
    if strcmp(varargin(1),{'plotname'}) == 1
        print('-dpdf',[varargin{2},num2str(plot_nr)])
    end
end

%%  Plot UKF Bias with sensor displacement and 1 RK substep

load(test_name)

s_b=1;
for ii=1:length(data(:,75))
    if data(ii,75)>0 && s_b==1
        eclip(1)=ii;
        s_b=0;
    end
    if data(ii,75)<1 && s_b==0
        eclip(2)=ii-1;
        break;
    end
end

stage(1)=max(max(abs(data(1:1000,50:52))));
stage(2)=max(max(abs(data(1000+1:2100,50:52))));
stage(3)=max(max(abs(data(2101:4095,50:52))));
stage(4)=max(max(abs(data(4096:6000,50:52))));
fprintf('\n  8%12.2f [deg]%12.2f [deg]%12.2f [deg]%12.2f [deg]\n',stage(1:4));

plot_nr=4;
fig_size=[6,21]; % cm
figure('PaperUnits', 'centimeters','PaperPosition', [0.0,(29.2-fig_size(1)) fig_size(2) fig_size(1)]); %210 � 297 A4
set(gcf,'DefaultAxesColorOrder',[0,0,0;0.6,0.6,0.6],'DefaultAxesLineStyleOrder','-|--');
set(gcf,'DefaultAxesFontSize',8,'DefaultAxesFontName', 'Times New Roman');
subplot(1,1,1);
plot(data(:,1),data(:,50:52),[eclip(1);eclip(2)],[-9 -9],'k-',[eclip(1);eclip(1)],[-9.5 -8.5],'k-',[eclip(2);eclip(2)],[-9.5 -8.5],'k-','LineWidth',0.9);
text(eclip(1)+100,-8,'Eclipse','FontSize',8,'FontName','Times New Roman');
ylabel('Error in Estimated Attitude [deg]');
xlabel('Time [s]');
legend('x','y','z')
axis([0 6000 -10 10]);
grid;

if nargin > 3
    if strcmp(varargin(1),{'plotname'}) == 1
        print('-dpdf',[varargin{2},num2str(plot_nr)])
    end
end


%% Check data
if check_data==1
    figure;
    plot(data(:,1),data(:,2:4));
    title('Real Attitude')
    grid;

    figure;
    plot(data(:,1),data(:,5:7));
    title('Real Angular Velocity')
    grid;

    figure;
    plot(data(:,1),data(:,8:10));
    title('Estimated Attitude (UKF)')
    grid;

    figure;
    plot(data(:,1),data(:,11:13));
    title('Estimated Angular Velocity (UKF)')
    grid;

    figure;
    plot(data(:,1),data(:,14:16));
    title('Residual Sun (UKF)')
    grid;

    figure;
    plot(data(:,1),data(:,17:19));
    title('Residual Mag (UKF)')
    grid;

    figure;
    plot(data(:,1),data(:,20:22));
    title('Residual Angular Velocity (UKF)')
    grid;

    figure;
    plot(data(:,1),data(:,23:25));
    title('Error in Attitude Estimation (UKF)')
    grid;

    figure;
    plot(data(:,1),data(:,26:28));
    title('Error in Estimated Angular Velocity (UKF)')
    grid;

    figure;
    plot(data(:,1),data(:,29:31));
    title('Estimated Attitude (UKF Bias)')
    grid;

    figure;
    plot(data(:,1),data(:,32:34));
    title('Estimated Angular Velocity (UKF Bias)')
    grid;

    figure;
    plot(data(:,1),data(:,35:37));
    title('Estimated Magnetometer Bias (UKF Bias)')
    grid;

    figure;
    plot(data(:,1),data(:,38:40));
    title('Estimated Sun Sensor Bias (UKF Bias)')
    grid;

    figure;
    plot(data(:,1),data(:,41:43));
    title('Residual Sun (UKF Bias)')
    grid;

    figure;
    plot(data(:,1),data(:,44:46));
    title('Residual Mag (UKF Bias)')
    grid;

    figure;
    plot(data(:,1),data(:,47:49));
    title('Residual Angular Velocity (UKF Bias)')
    grid;

    figure;
    plot(data(:,1),data(:,50:52));
    title('Error in Attitude Estimation (UKF Bias)')
    grid;

    figure;
    plot(data(:,1),data(:,53:55));
    title('Error in Estimated Angular Velocity (UKF Bias)')
    grid;

    figure;
    plot(data(:,1),data(:,56:58));
    title('Estimated Attitude (SVD)')
    grid;

    figure;
    plot(data(:,1),data(:,59:61));
    title('Error in Attitude Estimation (SVD)')
    grid;

    figure;
    plot(data(:,1),data(:,62));
    title('Sum of Covariance Matrix (SVD)')
    grid;

    figure;
    plot(data(:,1),data(:,63));
    title('Algortihm Execution Time (UKF)')
    grid;

    figure;
    plot(data(:,1),data(:,64));
    title('Algortihm Execution Time (UKF Bias)')
    grid;

    figure;
    plot(data(:,1),data(:,65));
    title('Algortihm Execution Time (SVD)')
    grid;

    figure;
    plot(data(:,1),data(:,66:68));
    title('Disturbance Torque')
    grid;

    figure;
    plot(data(:,1),data(:,69:71));
    title('Permanent Magnet Disturbance Torque')
    grid;

    figure;
    plot(data(:,1),data(:,72:74));
    title('Control Torque')
    grid;

    figure;
    plot(data(:,1),data(:,75));
    title('Eclipse Indicator (Real)')
    grid;

    figure;
    plot(data(:,1),data(:,76));
    title('Eclipse Indicator (SW)')
    grid;
end
end