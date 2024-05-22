% This file runs two simulations of the difference between the
% robust controlled NSO and the NSO controlled by simple state feedback. 
% Simulation params: q_init   = [-0.1166 -0.4432 0.4313 -0.7771] random
%                    q_ref    = [-0.1059 -0.7028 0.2078 -0.6721] generated
%                    omega_init = [0 0 0]
%                    drag_coef  = 2
%                    ro_air     = 9.25e-13 
%              sim1: I_sat      = [42.3 42.3 28.4]/1e3
%              sim2: I_sat      = 1.05*[42.3 42.3 28.4]/1e3

%clear all;
clc
hep = strrep(which('igrf2005.d','-all'),'igrf2005.d','');
cd(hep{1})

print_plots = true;
simulate = false;true;
outfile = 'main_test_1';

q_init   = [-0.1166 -0.4432 0.4313 -0.7771];
q_init   = q_init/norm(q_init);
q_ref    = [-0.1059 -0.7028 0.2078 -0.6721];
q_ref    = q_ref/norm(q_ref);
omega_init = [0 0 0];
drag_coef  = 2;
ro_air     = 9.25e-13;
ff_on      = 1;
t_init     = 100;
t_step     = 500;
t_sim      = 1000;

I_sat_unpert = [42.3 42.3 28.4]/1e3;
I_sat_pert = 0.95*[42.3 42.3 28.4]/1e3;
if simulate 
    tic
    %Simulation one - nominal Inertia matrix
    I_sat = I_sat_unpert;
    [T_main_all_1,X,Y_main_all_1] = sim('test/main_test_all',t_sim);
    [T_main_pp_1,X,Y_main_pp_1] = sim('test/main_test_pp',t_sim);
    
    %Simulation two - uncertainty + 5% on inertia matrix
    I_sat = I_sat_pert;
    [T_main_all_2,X,Y_main_all_2] = sim('test/main_test_all',t_sim);
    [T_main_pp_2,X,Y_main_pp_2] = sim('test/main_test_pp',t_sim);
    
    sim_time = toc;
    save(outfile, 'T_main*', 'Y_main*');
else
    load(outfile);
end
[bw_rob bw_sim] = calbw(diag(I_sat_unpert));
[bw_rob_pert bw_sim_pert] = calbw(diag(I_sat_pert));
disp('Bandwidth of norminal system:')
disp(strcat([9],'Robust system:',[32],[32],num2str(bw_rob(1)),'[rad/s]',[32],[32],'(',num2str(bw_rob(2)),'[Hz])'))
disp(strcat([9],'Simple system:',[32],[32],num2str(bw_sim(1)),'[rad/s]',[32],[32],'(',num2str(bw_sim(2)),'[Hz])'))
disp(' ')
disp('Bandwidth of uncertain system:')
disp(strcat([9],'Robust system:',[32],[32],num2str(bw_rob_pert(1)),'[rad/s]',[32],[32],'(',num2str(bw_rob_pert(2)),'[Hz])'))
disp(strcat([9],'Simple system:',[32],[32],num2str(bw_sim_pert(1)),'[rad/s]',[32],[32],'(',num2str(bw_sim_pert(2)),'[Hz])'))

%%%%%% Plot %%%%%%%%%%%%%%%%%
t_plot = 1:t_sim;
quart_zoom = 1.0;
vel_zoom = 0.001;%5e-9;
for k_main=1:2
    t =  evalin('base',strcat('T_main_all_',num2str(k_main)));
    y = evalin('base',strcat('Y_main_all_',num2str(k_main)));
    t_pp = evalin('base',strcat('T_main_pp_',num2str(k_main)));
    y_pp = evalin('base',strcat('Y_main_pp_',num2str(k_main)));
    att_pic_name = strcat('main_test_att_',num2str(k_main));
    vel_pic_name = strcat('main_test_vel_',num2str(k_main));
    
    figure(k_main*2-1)
    %quarternions
    subplot(3,1,1)
    plot(t,y(:,1:4))
    title('Attitude - Robust Control')
    ylabel('Attitude [.]')
    axis([t_plot(1) t_plot(end) -quart_zoom quart_zoom])
    grid on
    subplot(3,1,2)
    plot(t_pp,y_pp(:,1:4))
    title('Attitude - Simple Control')
    ylabel('Attitude [.]')
    axis([t_plot(1) t_plot(end) -quart_zoom quart_zoom])
    grid on
    subplot(3,1,3)
    plot(t,y(:,10),'b',t_pp,y_pp(:,10),'r')
    title('Attitude Error on z-axis')
    xlabel('Simulation time [s]')
    ylabel('Error angel [\circ]')
    axis([t_plot(1) t_plot(end) -0 1.6])
    legend('Robust','Simple')
    grid on
    if print_plots
        print(figure(k_main*2-1), '-depsc2', att_pic_name)
    end
    
    figure(k_main*2)
    %angular velocities
    subplot(2,1,1)
    plot(t,y(:,5:7))
    title('Angular Velocity - Robust Control')
    ylabel('Angular velocity [rad/s]')
    axis([t_plot(1) t_plot(end) -vel_zoom vel_zoom])
    grid on
    subplot(2,1,2)
    plot(t_pp,y_pp(:,5:7))
    title('Angular Velocity - Simple Control')
    xlabel('Simulation time [s]')
    ylabel('Angular velocity [rad/s]')
    axis([t_plot(1) t_plot(end) -vel_zoom vel_zoom])
    grid on
    if print_plots
        print(figure(k_main*2), '-depsc2', vel_pic_name)
    end
end

%Att in report
figure(11)
subplot(2,1,1)
t =  evalin('base',strcat('T_main_all_',num2str(1)));
y = evalin('base',strcat('Y_main_all_',num2str(1)));
t_pp = evalin('base',strcat('T_main_pp_',num2str(1)));
y_pp = evalin('base',strcat('Y_main_pp_',num2str(1)));
plot(t,y(:,10),'b',t_pp,y_pp(:,10),'r')
title('Attitude Error - Nominal Inertia')
ylabel('Error angel [\circ]')
legend('Robust','Simple')
axis([t_plot(1) t_plot(end) -0 1.6])
grid on
subplot(2,1,2)
t =  evalin('base',strcat('T_main_all_',num2str(2)));
y = evalin('base',strcat('Y_main_all_',num2str(2)));
t_pp = evalin('base',strcat('T_main_pp_',num2str(2)));
y_pp = evalin('base',strcat('Y_main_pp_',num2str(2)));
plot(t,y(:,10),'b',t_pp,y_pp(:,10),'r')
title('Attitude Error - Pertuped Inertia')
ylabel('Error angel [\circ]')
xlabel('Simulation time [s]')
axis([t_plot(1) t_plot(end) -0 1.6])
grid on
if print_plots
    print(figure(11), '-depsc2', 'main_test_1_angle')
end