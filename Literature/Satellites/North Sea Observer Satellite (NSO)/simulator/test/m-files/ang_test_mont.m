% This file is used to make a monte-carlo simulation of
% the angular velocity controller using different inertia
% matrices and random attitudes.
clear all
I_constant = [42.3 42.3 28.4]/1e3; % nominal inertia matrix
N_att      = 30;   % number of attitude changes for each inertia matrix
m_up_test  = 1.05; % upper limit of inertia matrix
m_do_test  = 0.95; % lower limit of inertia matrix

% defines the corners of the inertia matrix
corners_test = [m_do_test m_do_test m_do_test;
	  m_up_test m_do_test m_do_test;
	  m_up_test m_up_test m_do_test;
	  m_do_test m_up_test m_do_test;
	  m_do_test m_do_test m_up_test;
	  m_up_test m_do_test m_up_test;
	  m_up_test m_up_test m_up_test;
	  m_do_test m_up_test m_up_test;
	       1 1 1];

N = size(corners_test,1);

cdrho_min = 4e-14;
cdrho_max = 2*7.23e-12;
cdrho_delta = cdrho_max - cdrho_min;


sim_time = 0; % used to calculate the simulation time
tic 

file_name = 'ang_mont_estimation.mat'; %filename of mat file
momentum_dyn_on = 0; % momentum wheel dynamics off
est_on_off      = 1; % estimation and feedforward off (remember to change filename
ang_ref_value   = 0; % no reference is introduced
init_ang_vel    = [0 0 0]; % initial angular velocity
air_density     = 9.25e-13; % nominal air density
cd_value        = 1;
ref_change1     = 300; % reference step up time
ref_change2     = 2500;% reference step down time
hep_hey         = 0;   % count variable
for l_ang_mont=1:N
  inertial_matrix = diag(I_constant)*corners_test(l_ang_mont,:)';
  for k_ang_mont=1:N_att
    theta = rand*2*pi;
    phi = rand*2*pi;
    alpha = rand*2*pi;
    %Random e-vector
    [e(1) e(2) e(3)]=sph2cart(theta,phi,1);
    %Random quaternion used in the mdl mask of the NSO.
    %Random constant of Cd*rho
    air_density=cdrho_min+rand*cdrho_delta;
    ang_attitude = [e(1)*sin(alpha/2) e(2)*sin(alpha/2) e(3)*sin(alpha/2) ...
		    cos(alpha/2)];
    hep_hey = hep_hey + 1;
    [t_mont{hep_hey},x,y_mont{hep_hey}] = sim('test/angular_control.mdl');
    toc - sim_time
    sim_time = toc;
  end
end

save(file_name)








