% This file is used to simulated an orbit with the push-broom
% controller. The NSO is commanded to point towards Aalborg and
% the error between the reference quaternion and the actual
% attitude is plotted along with the attitude and angular
% velocity.

%clear all
%clc

%Clever nifty little thing that changes the dir to the one with igrf2005.d,
%such that one does not have to do it manually to run the sims
hep = strrep(which('igrf2005.d','-all'),'igrf2005.d','');
cd(hep{1})

tic
simulate = false;
plot_figures = true;
print_plots = false;

%Static sim parameters
omega_init = [0 0 0];
drag_coef  = 2;
ff_on      = 1;

operating_point_constants; %Someone clears I
I_sat = I*[1 1 1]';
q_init = [sqrt(0.5) 0 0 sqrt(0.5)];

ro_air = 9.25e-13;

ref_lat  = 9.54;
ref_long = 57.02;

R_earth  = 6.371e6;

[x,y,z]  = sph2cart(deg2rad(ref_lat),deg2rad(ref_long),R_earth);
R_ref_E =[x,y,z];

if simulate
  %Simulation on NSO robust configuration
  [T_vec_ref,X,Y_vec_ref]=sim('test/vector_reference.mdl');
end

sim_time = toc;
%%%%%%%%%%%%%%%% Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error_angle = Y_vec_ref(:,10);
omega_vec_ref = Y_vec_ref(:,5:7);
quat_vec_ref = Y_vec_ref(:,1:4);

if plot_figures
  
  figure(1)
  subplot(2,1,1)
  plot(T_vec_ref,quat_vec_ref)
  title('Tracking of the North Sea')
  xlabel('Simulation time [s]')
  ylabel('Quaternion')
  axis([0 T_vec_ref(end) -1 1])
  grid on
  
  subplot(2,1,2)
  plot(T_vec_ref,omega_vec_ref)
 % title('Angular Velocity')
  xlabel('Simulation time [s]')
  ylabel('Angular velocity [rad/s]')
  legend('x-axis','y-axis','z-axis')
  axis([0 T_vec_ref(end) -3e-3 3e-3])
  grid on

  if print_plots
    print -depsc2 'vector_reference_w_q.eps'
  end
  figure(2)
  plot(T_vec_ref,error_angle)
  title('Attitude Error on z-axis')
  xlabel('Simulation time [s]')
  ylabel('Error angle [\circ]')
  axis([0 T_vec_ref(end) 0 3])
  grid on
  
  if print_plots
    print -depsc2 'vector_reference_z_error.eps'
  end
end
