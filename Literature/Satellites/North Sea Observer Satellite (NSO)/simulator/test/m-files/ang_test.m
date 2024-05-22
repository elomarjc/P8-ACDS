% This file is used to generate the plots used in the
% evaluation of the angular velocity controller. It used the
% model file angular_control.mdl which is located in the
% test folder in simulation.

simulate_sim = true;
print_plots = 1;   % determines if the plots should be saved
ref_change1 = 300; % time of step in reference
ref_change2 = 2500;% time of 2nd step in reference


ang_control_config = 'ang_robust_control.m';
ang_estimation_config = 'kalman_estimator.m';
if simulate_sim
  tic
  % simulation with momentum wheel dynamics
  % estimation and feedforward OFF
  momentum_dyn_on = 1; % momentum wheel dynamics on
  est_on_off      = 0; % estimation and feedforward off
  ang_ref_value   = 0; % no reference is introduced
  init_ang_vel    = [0.0021 -0.0021 0.0021];
  ang_attitude    = [sqrt(0.5) 0 0 sqrt(0.5)];
  inertial_matrix = [42.3 42.3 28.4]/1e3;
  air_density     = 9.25e-13;
  cd_value        = 2;
  [t_est_off,x,y_est_off] = sim('test/angular_control.mdl');
  
  toc

  % Add a breakpoint here
  % keyboard;
  
  % simulation with momentum wheel dynamics included
  % estimation and feedforward on
  momentum_dyn_on = 1; % momentum wheel dynamics on
  est_on_off      = 1; % estimation and feedforward on
  ang_ref_value   = 0; % no reference is introduced
  init_ang_vel    = [0.0021 -0.0021 0.0021];
  ang_attitude    = [sqrt(0.5) 0 0 sqrt(0.5)];
  inertial_matrix = [42.3 42.3 28.4]/1e3;
  air_density     = 9.25e-13;
  cd_value        = 2;
  [t_est_on,x,y_est_on] = sim('test/angular_control.mdl');
  
  toc

  % Add a breakpoint here
  % keyboard;
  
  % simulation with momentum wheel dynamics included
  % estimation and feedforward on
  % reference input
  momentum_dyn_on = 1; % momentum wheel dynamics on
  est_on_off      = 1; % estimation and feedforward on
  ang_ref_value   = 0.0045; % reference signal is 0.0045
  init_ang_vel    = [0.0021 -0.0021 0.0021];
  ang_attitude    = [sqrt(0.5) 0 0 sqrt(0.5)];
  inertial_matrix = [42.3 42.3 28.4]/1e3;
  air_density     = 9.25e-13;
  cd_value        = 2;
  [t_ref,x,y_ref] = sim('test/angular_control.mdl');
  
  toc

  % Add a breakpoint here
  % keyboard;
  
  % simulation withOUT momentum wheel dynamics included
  % estimation and feedforward on
  momentum_dyn_on = 0; % momentum wheel dynamics off
  est_on_off      = 1; % estimation and feedforward on
  ang_ref_value   = 0; % no reference is introduced
  init_ang_vel    = [0.0021 -0.0021 0.0021];
  ang_attitude    = [sqrt(0.5) 0 0 sqrt(0.5)];
  inertial_matrix = [42.3 42.3 28.4]/1e3;
  air_density     = 9.25e-13;
  cd_value        = 2;
  [t_no_dyn,x,y_no_dyn] = sim('test/angular_control.mdl');
  
  toc

  % Add a breakpoint here
  % keyboard;
  
  
  % simulation with momentum wheel dynamics included
  % estimation and feedforward on - simple system
  ang_control_config = 'basic_ang_control.m';
  ang_estimation_config = 'linear_estimator.m';
  
  momentum_dyn_on = 1; % momentum wheel dynamics on
  est_on_off      = 1; % estimation and feedforward on
  ang_ref_value   = 0; % no reference is introduced
  init_ang_vel    = [0.0021 -0.0021 0.0021];
  ang_attitude    = [sqrt(0.5) 0 0 sqrt(0.5)];
  inertial_matrix = [42.3 42.3 28.4]/1e3;
  air_density     = 9.25e-13;
  cd_value        = 2;
  [t_lin,x,y_lin] = sim('test/angular_control.mdl');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots the simulation results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
%estimation off
plot(t_est_off,y_est_off(:,1:3))
title('Angular Velocity Controller without Feedforward')
xlabel('Simulation time [s]')
ylabel('Angular velocity [rad/s]')
legend('x-axis','y-axis','z-axis')
xlim([0 t_est_off(end)])
grid on
axis([t_est_off(1) t_est_off(end) -1.5e-6 1.5e-6])

if print_plots
    print -depsc2 'ang_control_est_off.png'
end

figure(2)
%estimation off - zoom
hephey = find(t_est_off < 20);
plot(t_est_off(hephey),y_est_off(hephey,1:3))
title('Angular Velocity Controller without Feedforward - Settling')
xlabel('Simulation time [s]')
ylabel('Angular velocity [rad/s]')
legend('x-axis','y-axis','z-axis')
xlim([0 t_est_off(hephey(end))])
grid on
%axis([t_est_off(find(t_est_off < 20)) -2e-3 5e-3])

if print_plots
    print -depsc2 'ang_control_est_off_zoom.png'
end

figure(3)
%estimation on
plot(t_est_on,y_est_on(:,1:3))
title('Angular Velocity Controller with Feedforward')
xlabel('Simulation time [s]')
ylabel('Angular velocity [rad/s]')
legend('x-axis','y-axis','z-axis')
grid on
axis([t_est_on(1) t_est_on(end) -3e-7 3e-7])

if print_plots
    print -depsc2 'ang_control_est_on.png'
end

figure(4)
%estimation on - zoom
hephey = find(t_est_on < 20);
plot(t_est_on(hephey),y_est_on(hephey,1:3))
title('Angular Velocity Controller with Feedforward - Settling')
xlabel('Simulation time [s]')
ylabel('Angular velocity [rad/s]')
legend('x-axis','y-axis','z-axis')
xlim([0 t_est_on(hephey(end))])
grid on

if print_plots
    print -depsc2 'ang_control_est_on_zoom.eps'
end

figure(5)
%momentum dyn on
hephey = find(t_no_dyn < 20);
plot(t_no_dyn(hephey),y_no_dyn(hephey,1:3))
title('Angular Velocity Controller without Actuator Dynamics')
xlabel('Simulation time [s]')
ylabel('Angular velocity [rad/s]')
legend('x-axis','y-axis','z-axis')
grid on
xlim([0 t_no_dyn(hephey(end))])

if print_plots
    print -depsc2 'ang_control_no_dyn.eps'
end

figure(6)
%estimation on Simple
plot(t_lin,y_lin(:,1:3))
title('Simple Angular Velocity Controller with Feedforward')
xlabel('Simulation time [s]')
ylabel('Angular velocity [rad/s]')
legend('x-axis','y-axis','z-axis')
grid on
axis([t_lin(1) t_lin(end) -3e-7 3e-7])

if print_plots
    print -depsc2 'ang_control_est_on_simple.eps'
end

figure(7)
%with reference input
subplot(2,1,1)
plot(t_ref,y_ref(:,1:3))
title('Angular Velocity Vontroller with Reference Input')
xlabel('Simulation time [s]')
ylabel('Angular velocity [rad/s]')
legend('x-axis','y-axis','z-axis')
grid on
axis([t_ref(1) t_ref(end) -1e-3 5e-3])
subplot(2,2,3)
hephey = find(295 <t_ref & t_ref < 315);
plot(t_ref(hephey),y_ref(hephey,1:3))
xlabel('Simulation time [s]')
ylabel('Angular velocity [rad/s]')
grid on
xlim([t_ref(hephey(1)) t_ref(hephey(end))])

subplot(2,2,4)
hephey = find(2495 <t_ref & t_ref < 2515);
plot(t_ref(hephey),y_ref(hephey,1:3))
xlabel('Simulation time [s]')
ylabel('Angular velocity [rad/s]')
grid on
xlim([t_ref(hephey(1)) t_ref(hephey(end))])


if print_plots
    print -depsc2 'ang_control_reference.eps'
end

