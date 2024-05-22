% This file runs a monte carlo simulation of the difference between the
% estimated disturbance and the actual disturbance. The parameters that
% vary are the inertia matrix, Cd and rho. Every call of the mdl file runs
% a simulation specified to take one orbit = 6179s.
tic
sim_time = 0;
simulate = false;
do_est_err_plot = true;
print_plots = true;
%Variation of inertia matrix when deployed.
Ivar = 5; %variation in percent
I1d = 42.3/1e3;
I2d = 42.3/1e3;
I3d = 28.4/1e3;

I1dmin = I1d - (I1d/100)*Ivar;
I2dmin = I2d - (I2d/100)*Ivar;
I3dmin = I3d - (I3d/100)*Ivar;

I1dmax = I1d + (I1d/100)*Ivar;
I2dmax = I2d + (I2d/100)*Ivar;
I3dmax = I3d + (I3d/100)*Ivar;

I1delta = I1dmax-I1dmin;
I2delta = I2dmax-I2dmin;
I3delta = I3dmax-I3dmin;

%Inertia matrices for nominal situation and the 8 vertices.
Inom = [I1d I2d I3d];
Iver1 = [I1dmin I2dmin I3dmin];
Iver2 = [I1dmax I2dmin I3dmin];
Iver3 = [I1dmax I2dmax I3dmin];
Iver4 = [I1dmin I2dmax I3dmin];
Iver5 = [I1dmin I2dmin I3dmax];
Iver6 = [I1dmax I2dmin I3dmax];
Iver7 = [I1dmax I2dmax I3dmax];
Iver8 = [I1dmin I2dmax I3dmax];

I = [Inom Iver1 Iver2 Iver3 Iver4 Iver5 Iver6 Iver7 Iver8];

%Variation of Cd*rho
cdrho_min = 4e-14;
cdrho_max = 2*7.23e-12;
cdrho_delta = cdrho_max - cdrho_min;
if ~do_est_err_plot
  dist_N = 39;
  dist_k = 30;
else
  dist_N = 1;
  dist_k = 1;
end
if simulate
  for N=1:dist_N
    %Generating random inertia matrix, Id is included in the mdl mask of
    %the NSO.
    if(N>30)
      Id = I(N-30);
    else
      Id=[(I1dmin+rand*I1delta) (I2dmin+rand*I2delta) (I3dmin+rand*I3delta)];
    end
    for k=1:dist_k
      theta = rand*2*pi;
      phi = rand*2*pi;
      alpha = rand*2*pi;
      %Random e-vector
      [e(1) e(2) e(3)]=sph2cart(theta,phi,1);
      %Random quaternion used in the mdl mask of the NSO.
      q_att = [e(1)*sin(alpha/2) e(2)*sin(alpha/2) e(3)*sin(alpha/2) cos(alpha/2)];
      %Random constant of Cd*rho, cdrho is used in the mdl mask of the
      %NSO.
      cdrho=cdrho_min+rand*cdrho_delta;
      %Running 30 simulations with one inertia matrix and 30 different
      %starting attitudes and constants.
      [T,X,Y{k}]=sim('test/uncertainty_disturbance_torques');
      toc - sim_time
      sim_time = toc;
    end
    sim_dist{N}=Y;
  end
end
if ~do_est_err_plot & simulate
  save dist_sim sim_dist;        
else
  figure(1)
  plot(T,Y{1})
  xlabel('Simulation time [s]')
  ylabel('Disturbance residual [Nm]')
  title('Estimation Error')
  legend('x-axis','y-axis','z-axis')
  xlim([0 T(end)])
  grid on
  %axis([0 T(end) min(Y{1}) max(Y{1})])
  if print_plots
    print -depsc2 'ang_control_estimation_error.eps'
  end
end