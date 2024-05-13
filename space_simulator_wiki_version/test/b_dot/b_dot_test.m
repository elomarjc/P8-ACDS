clear
close all
clc

Simulate=0; % Run simulations 1=true,0=false
Plot_sim=1; % Make plots 1=true,0=false
Extra_stuff=0; % Calculate extra stuff 1=true,0=false

save_dir='test/b_dot/';

%% Angular difference between nadir vector and B-field vector over Aalborg
if Extra_stuff==1
    % B-field vector over aalborg
    % http://www.ngdc.noaa.gov/geomagmodels/struts/calcIGRFWMM
    V_b=[12998.5,45.3,36231.7];
    V_b_norm=V_b/norm(V_b);
    
    % Nadir vector
    V_n=[0,0,1];
    
    % Angle in degree
    angle_z_B_aalborg=180/pi*acos(V_b_norm*V_n')
end

%% Bodeplots of filter
if Extra_stuff==1
    s=tf('s');
    wc=0.7;
    H=wc*s/(s+wc);
    figure
    bode(H)
    G=c2d(H,0.1);
    figure
    bode(G)
    
    s=tf('s');
    wc=0.7*10000;
    H=wc*s/(s+wc);
    figure
    bode(H)
    G=c2d(H,0.1);
    figure
    bode(G)
end

%% Filter simulations (angular_vel=10)
if Extra_stuff==1
    % variables
    sim_duration=200;
    angular_vel=10;
    temp_coil=85;

    sim('b_dot_filter', sim_duration);
    
    fig_size=[12,18]; % cm
    figure('PaperUnits', 'centimeters','PaperPosition', [0.0,(29.2-fig_size(1)) fig_size(2) fig_size(1)]); %210 � 297 A4
    set(gcf,'DefaultAxesFontSize',8,'DefaultAxesFontName', 'Times New Roman');
    set(gcf,'DefaultAxesColorOrder',[0,0,0;0.6,0.6,0.6],'DefaultAxesLineStyleOrder','-|--');
    
    subplot(2,2,1)
    plot(sim_time,bdot_c_matlab)
    ylabel('Continous B-dot (dx/dt simulink blok) [T/s]')
    legend('Bdot_x','Bdot_y','Bdot_z')
    axis([0 sim_duration -8*10^(-6) 8*10^(-6)]);
    grid;
    
    subplot(2,2,2)
    plot(sim_time,bdot_c_filter)
    ylabel('Continous B-dot (filter) [T/s]')
    legend('Bdot_x','Bdot_y','Bdot_z')
    axis([0 sim_duration -8*10^(-6) 8*10^(-6)]);
    grid;
    
    subplot(2,2,3)
    h=get(gca,'Position');
    set(gca,'Position',[h(1) h(2)+0.04 h(3) h(4)])
    time_disc=[0:0.1:sim_duration];
    plot(time_disc,bdot_d_nofilter)
    ylabel('Discrete B-dot (no filter) [T/s]')
    xlabel('Time [s]')
    legend('Bdot_x','Bdot_y','Bdot_z')
    axis([0 sim_duration -8*10^(-5) 8*10^(-5)]);
    grid;
    
    subplot(2,2,4)
    h=get(gca,'Position');
    set(gca,'Position',[h(1) h(2)+0.04 h(3) h(4)])
    plot(time_disc,bdot_d_filter)
    ylabel('Discrete B-dot (filter) [T/s]')
    xlabel('Time [s]')
    legend('Bdot_x','Bdot_y','Bdot_z')
    axis([0 sim_duration -8*10^(-6) 8*10^(-6)]);
    grid; 
    print('-dpdf','bdot_filter_sim')
end

%% B-dot simulation (temp_coil=85,angular_vel=10)
if Simulate==1
    fprintf('B-dot simulation 1 \n');

    % variables
    sim_duration=4*5830;
    angular_vel=10;
    temp_coil=85;
    w_c=0.7;
    
    sim('b_dot_sim', sim_duration); 
    
    save([save_dir,'b_dot_temp85_omega10_data','.mat'],'temp_coil','angular_vel','sim_duration','w_c','omega_s','power_mt','angle_z_B');
end

%% B-dot simulation (temp_coil=30,angular_vel=10)
if Simulate==1
    fprintf('B-dot simulation 2 \n');

    % variables
    sim_duration=4*5830;
    angular_vel=10;
    temp_coil=30;
    w_c=0.7;
    
    sim('b_dot_sim', sim_duration); 
    
    save([save_dir,'b_dot_temp30_omega10_data','.mat'],'temp_coil','angular_vel','sim_duration','w_c','omega_s','power_mt','angle_z_B');
end

%% B-dot simulation (temp_coil=-25,angular_vel=10)
if Simulate==1
    fprintf('B-dot simulation 3 \n');

    % variables
    sim_duration=4*5830;
    angular_vel=10;
    temp_coil=-25;
    w_c=0.7;
    
    sim('b_dot_sim', sim_duration); 
    
    save([save_dir,'b_dot_tempm25_omega10_data','.mat'],'temp_coil','angular_vel','sim_duration','w_c','omega_s','power_mt','angle_z_B');
end

%% B-dot simulation (temp_coil=30,angular_vel=500)
if Simulate==1
    fprintf('B-dot simulation 4 \n');
    
    % variables
    sim_duration=8*5830;
    angular_vel=500;
    temp_coil=30;
    w_c=0.7*10000;
    
    sim('b_dot_sim', sim_duration); 
    
    save([save_dir,'b_dot_temp30_omega500_data','.mat'],'temp_coil','angular_vel','sim_duration','w_c','omega_s','power_mt','angle_z_B');
end

%% Plot figures
if Plot_sim==1
    clear
    Dutycycle=0.88;
    Ts_bdot=0.1;
    
    load b_dot_temp85_omega10_data
    
    % Detumble time in [orbits]
    for i=1:sim_duration
        if omega_s(i,1)>0.3 || omega_s(i,1)<-0.3 || omega_s(i,2)>0.3 || omega_s(i,2)<-0.3 || omega_s(i,3)>0.3 || omega_s(i,3)<-0.3
            detumble_time_sim1=i;
        end
    end
    detumble_time_sim1=detumble_time_sim1/5830
    
    % power calc in [W]
    Time_start=0;
    Time_end=5830;
    Avrg_Power_per_axis_sim1_orbit1=sum(power_mt(Time_start/Ts_bdot+1:Time_end/Ts_bdot+1,:))*Dutycycle/((Time_end-Time_start)/Ts_bdot);
    Avrg_Power_total_sim1_orbit1=sum(Avrg_Power_per_axis_sim1_orbit1)+6*0.005*Dutycycle
    Time_start=5830;
    Time_end=2*5830;
    Avrg_Power_per_axis_sim1_orbit2=sum(power_mt(Time_start/Ts_bdot+1:Time_end/Ts_bdot+1,:))*Dutycycle/((Time_end-Time_start)/Ts_bdot);
    Avrg_Power_total_sim1_orbit2=sum(Avrg_Power_per_axis_sim1_orbit2)+6*0.005*Dutycycle
    Time_start=2*5830;
    Time_end=3*5830;
    Avrg_Power_per_axis_sim1_orbit3=sum(power_mt(Time_start/Ts_bdot+1:Time_end/Ts_bdot+1,:))*Dutycycle/((Time_end-Time_start)/Ts_bdot);
    Avrg_Power_total_sim1_orbit3=sum(Avrg_Power_per_axis_sim1_orbit3)+6*0.005*Dutycycle
    Time_start=3*5830;
    Time_end=4*5830;
    Avrg_Power_per_axis_sim1_orbit4=sum(power_mt(Time_start/Ts_bdot+1:Time_end/Ts_bdot+1,:))*Dutycycle/((Time_end-Time_start)/Ts_bdot);
    Avrg_Power_total_sim1_orbit4=sum(Avrg_Power_per_axis_sim1_orbit4)+6*0.005*Dutycycle
    
    load b_dot_temp30_omega10_data
    
    % Detumble time in [orbits]
    for i=1:sim_duration
        if omega_s(i,1)>0.3 || omega_s(i,1)<-0.3 || omega_s(i,2)>0.3 || omega_s(i,2)<-0.3 || omega_s(i,3)>0.3 || omega_s(i,3)<-0.3
            detumble_time_sim2=i;
        end
    end
    detumble_time_sim2=detumble_time_sim2/5830
    
    fig_size=[18,18]; % cm
    figure('PaperUnits', 'centimeters','PaperPosition', [0.0,(29.2-fig_size(1)) fig_size(2) fig_size(1)]); %210 � 297 A4
    set(gcf,'DefaultAxesFontSize',8,'DefaultAxesFontName', 'Times New Roman');
    set(gcf,'DefaultAxesColorOrder',[0,0,0;0.6,0.6,0.6],'DefaultAxesLineStyleOrder','-|--');
    subplot(3,1,1)
    time_plot=0:((sim_duration/5830)/sim_duration):(sim_duration/5830);
    plot(time_plot,omega_s)
    set(gca,'ytick',(-15:5:15))
    ylabel('Angular Velocity [deg/s]')
    legend('^s\omega_x','^s\omega_y','^s\omega_z','Location','SouthEast')
    axis([0 sim_duration/5830 -15 15]);
    grid;
    subplot(3,1,2)
    h=get(gca,'Position');
    set(gca,'Position',[h(1) h(2)+0.04 h(3) h(4)])
    time_plot2=2:((sim_duration/5830)/sim_duration):3;
    plot(time_plot2,omega_s(2*5830:3*5830,:))
    set(gca,'ytick',(-0.3:0.1:0.3))
    ylabel('Angular Velocity [deg/s]')
    legend('^s\omega_x','^s\omega_y','^s\omega_z','Location','SouthEast')
    axis([2 3 -0.3 0.3]);
    grid;
    subplot(3,1,3)
    h=get(gca,'Position');
    set(gca,'Position',[h(1) h(2)+0.08 h(3) h(4)])
    plot(time_plot,angle_z_B)
    set(gca,'ytick',(0:30:180))
    ylabel('Angle [deg]')
    xlabel('Time [orbits]')
    axis([0 sim_duration/5830 0 180]);
    legend('Angular deviation between SBRF z-axis and B-field vector','Location','SouthEast')
    grid;
    print('-dpdf','bdot_plot1')
    
    % power calc in [W]
    Time_start=0;
    Time_end=5830;
    Avrg_Power_per_axis_sim2_orbit1=sum(power_mt(Time_start/Ts_bdot+1:Time_end/Ts_bdot+1,:))*Dutycycle/((Time_end-Time_start)/Ts_bdot);
    Avrg_Power_total_sim2_orbit1=sum(Avrg_Power_per_axis_sim2_orbit1)+6*0.005*Dutycycle
    Time_start=5830;
    Time_end=2*5830;
    Avrg_Power_per_axis_sim2_orbit2=sum(power_mt(Time_start/Ts_bdot+1:Time_end/Ts_bdot+1,:))*Dutycycle/((Time_end-Time_start)/Ts_bdot);
    Avrg_Power_total_sim2_orbit2=sum(Avrg_Power_per_axis_sim2_orbit2)+6*0.005*Dutycycle
    Time_start=2*5830;
    Time_end=3*5830;
    Avrg_Power_per_axis_sim2_orbit3=sum(power_mt(Time_start/Ts_bdot+1:Time_end/Ts_bdot+1,:))*Dutycycle/((Time_end-Time_start)/Ts_bdot);
    Avrg_Power_total_sim2_orbit3=sum(Avrg_Power_per_axis_sim2_orbit3)+6*0.005*Dutycycle
    Time_start=3*5830;
    Time_end=4*5830;
    Avrg_Power_per_axis_sim2_orbit4=sum(power_mt(Time_start/Ts_bdot+1:Time_end/Ts_bdot+1,:))*Dutycycle/((Time_end-Time_start)/Ts_bdot);
    Avrg_Power_total_sim2_orbit4=sum(Avrg_Power_per_axis_sim2_orbit4)+6*0.005*Dutycycle
    
    load b_dot_tempm25_omega10_data
    
    % Detumble time in [orbits]
    for i=1:sim_duration
        if omega_s(i,1)>0.3 || omega_s(i,1)<-0.3 || omega_s(i,2)>0.3 || omega_s(i,2)<-0.3 || omega_s(i,3)>0.3 || omega_s(i,3)<-0.3
            detumble_time_sim3=i;
        end
    end
    detumble_time_sim3=detumble_time_sim3/5830
    
    % power calc in [W]
    Time_start=0;
    Time_end=5830;
    Avrg_Power_per_axis_sim3_orbit1=sum(power_mt(Time_start/Ts_bdot+1:Time_end/Ts_bdot+1,:))*Dutycycle/((Time_end-Time_start)/Ts_bdot);
    Avrg_Power_total_sim3_orbit1=sum(Avrg_Power_per_axis_sim3_orbit1)+6*0.005*Dutycycle
    Time_start=5830;
    Time_end=2*5830;
    Avrg_Power_per_axis_sim3_orbit2=sum(power_mt(Time_start/Ts_bdot+1:Time_end/Ts_bdot+1,:))*Dutycycle/((Time_end-Time_start)/Ts_bdot);
    Avrg_Power_total_sim3_orbit2=sum(Avrg_Power_per_axis_sim3_orbit2)+6*0.005*Dutycycle
    Time_start=2*5830;
    Time_end=3*5830;
    Avrg_Power_per_axis_sim3_orbit3=sum(power_mt(Time_start/Ts_bdot+1:Time_end/Ts_bdot+1,:))*Dutycycle/((Time_end-Time_start)/Ts_bdot);
    Avrg_Power_total_sim3_orbit3=sum(Avrg_Power_per_axis_sim3_orbit3)+6*0.005*Dutycycle
    Time_start=3*5830;
    Time_end=4*5830;
    Avrg_Power_per_axis_sim3_orbit4=sum(power_mt(Time_start/Ts_bdot+1:Time_end/Ts_bdot+1,:))*Dutycycle/((Time_end-Time_start)/Ts_bdot);
    Avrg_Power_total_sim3_orbit4=sum(Avrg_Power_per_axis_sim3_orbit4)+6*0.005*Dutycycle
    
    load b_dot_temp30_omega500_data
    
    % Detumble time in [orbits]
    for i=1:sim_duration
        if omega_s(i,1)>0.3 || omega_s(i,1)<-0.3 || omega_s(i,2)>0.3 || omega_s(i,2)<-0.3 || omega_s(i,3)>0.3 || omega_s(i,3)<-0.3
            detumble_time_sim4=i;
        end
    end
    detumble_time_sim4=detumble_time_sim4/5830
    
    fig_size=[6,18]; % cm
    figure('PaperUnits', 'centimeters','PaperPosition', [0.0,(29.2-fig_size(1)) fig_size(2) fig_size(1)]); %210 � 297 A4
    set(gcf,'DefaultAxesFontSize',8,'DefaultAxesFontName', 'Times New Roman');
    set(gcf,'DefaultAxesColorOrder',[0,0,0;0.6,0.6,0.6],'DefaultAxesLineStyleOrder','-|--');
    time_plot=0:((sim_duration/5830)/sim_duration):(sim_duration/5830);
    plot(time_plot,omega_s)
    set(gca,'ytick',(-800:200:800))
    ylabel('Angular Velocity [deg/s]')
    xlabel('Time [orbits]')
    legend('^s\omega_x','^s\omega_y','^s\omega_z','Location','SouthEast')
    axis([0 sim_duration/5830 -800 800]);
    grid;
    print('-dpdf','bdot_plot2')

    % power calc in [W]
    Time_start=0;
    Time_end=5830;
    Avrg_Power_per_axis_sim4_orbit1=sum(power_mt(Time_start/Ts_bdot+1:Time_end/Ts_bdot+1,:))*Dutycycle/((Time_end-Time_start)/Ts_bdot);
    Avrg_Power_total_sim4_orbit1=sum(Avrg_Power_per_axis_sim4_orbit1)+6*0.005*Dutycycle
    Time_start=5830;
    Time_end=2*5830;
    Avrg_Power_per_axis_sim4_orbit2=sum(power_mt(Time_start/Ts_bdot+1:Time_end/Ts_bdot+1,:))*Dutycycle/((Time_end-Time_start)/Ts_bdot);
    Avrg_Power_total_sim4_orbit2=sum(Avrg_Power_per_axis_sim4_orbit2)+6*0.005*Dutycycle
    Time_start=2*5830;
    Time_end=3*5830;
    Avrg_Power_per_axis_sim4_orbit3=sum(power_mt(Time_start/Ts_bdot+1:Time_end/Ts_bdot+1,:))*Dutycycle/((Time_end-Time_start)/Ts_bdot);
    Avrg_Power_total_sim4_orbit3=sum(Avrg_Power_per_axis_sim4_orbit3)+6*0.005*Dutycycle
    Time_start=3*5830;
    Time_end=4*5830;
    Avrg_Power_per_axis_sim4_orbit4=sum(power_mt(Time_start/Ts_bdot+1:Time_end/Ts_bdot+1,:))*Dutycycle/((Time_end-Time_start)/Ts_bdot);
    Avrg_Power_total_sim4_orbit4=sum(Avrg_Power_per_axis_sim4_orbit4)+6*0.005*Dutycycle
    
end