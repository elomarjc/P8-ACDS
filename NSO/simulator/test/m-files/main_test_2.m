% This file runs a monte carlo simulation of the difference between the
% robust controlled NSO and the NSO controlled by simple state feedback.
% The parameters that vary are the inertia matrix, Cd and rho (one var)
% and the initial attitude. Every call of the mdl file runs a one 
%orbit simulation.

%clear all
clc
%Clever nifty little thing that changes the dir to the one with igrf2005.d,
%such that one does not have to do it manually to run the sims
hep = strrep(which('igrf2005.d','-all'),'igrf2005.d','');
cd(hep{1})

tic
simulate = false;%true;
playsound = false;
print_plot = false;
%Static sim parameters
q_ref    = [0 0 0 1]; %Not used
omega_init = [0 0 0];
drag_coef  = 1;
ff_on      = 1;
t_step     = 0;
t_init     = 0;

%Variation of inertia matrix when deployed.
m_up = 1.05;
m_do = 0.95;

% vertices of the convex uncertainty region
uncert_reg = [m_do m_do m_do;
    m_up m_do m_do;
    m_up m_up m_do;
    m_do m_up m_do;
    m_do m_do m_up;
    m_up m_do m_up;
    m_up m_up m_up;
    m_do m_up m_up;
    1 1 1];

%Get inertia
operating_point_constants;
I_norm = I;

%Variation of Cd*rho - set Cd = 1
cdrho_min = 4e-14;
cdrho_max = 2*7.23e-12;
cdrho_delta = cdrho_max - cdrho_min;

N_main = size(uncert_reg,1);
%N_main = 2;
Q_main = 30;% Number of iterations with dif. att.


%Estimator settle time
settle_time = 150;
if simulate
    %Dummy initiation vars
    z_angle_all = [];
    z_angle_pp  = [];
    omega_all = [];
    omega_pp = [];
    for s_main=1:N_main
        % Nine inertia matrices (8 vertirces and 1 nominal)
        I_sat = I_norm*uncert_reg(s_main,:)';
        for k_main=1:Q_main
            disp(strcat('*** Iteration :',[32],num2str(s_main),' and :',[32],num2str(k_main))) 
            %Random e-vector
            theta = rand*2*pi;
            phi = rand*2*pi;
            [e(1) e(2) e(3)]=sph2cart(theta,phi,1);
            %Random quaternion used in the mdl mask of the NSO.
            alpha = rand*2*pi;
            q_init = [e(1)*sin(alpha/2) e(2)*sin(alpha/2) e(3)*sin(alpha/2) cos(alpha/2)];
            %Random constant of Cd*rho
            ro_air=cdrho_min+rand*cdrho_delta;
            %Simulation on NSO robust configuration
            [T_all{k_main},X,Y_all{k_main}]=sim('test/main_test_all');
            %Simulation on simple control of NSO
            [T_pp{k_main},X,Y_pp{k_main}]=sim('test/main_test_pp',T_all{k_main});
            %dont use settle time
            asdf = find(T_all{k_main} >= settle_time,1,'first');
            %save data
            z_angle_all = [z_angle_all Y_all{k_main}(asdf:end,10)'];
            z_angle_pp  = [z_angle_pp Y_pp{k_main}(asdf:end,10)'];
            omega_all = [omega_all Y_all{k_main}(asdf:end,5:7)'];
            omega_pp = [omega_pp Y_pp{k_main}(asdf:end,5:7)'];
        end
        time{s_main}=T_all;
        inertia_var_all{s_main}=Y_all;
        inertia_var_pp{s_main}=Y_pp;
    end
    save monte_main_test time* inertia_var* z_angle* omega*;
    sim_time = toc;
    if playsound
        [x,fs] = wavread('tada.wav');
        sound(x,fs);
    end
end

%%%%%%%%%%%%%%%% Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%statistics
z_all_mean = mean(rad2deg(z_angle_all));
z_all_sigma = std(rad2deg(z_angle_all));
z_all_max = max(rad2deg(z_angle_all));
z_all_min = min(rad2deg(z_angle_all));
omega_all_mean = mean(omega_all,2);
omega_all_sigma = std(omega_all,0,2);
omega_all_max = max(omega_all,[],2);
omega_all_min = min(omega_all,[],2);


z_pp_mean = mean(rad2deg(z_angle_pp));
z_pp_sigma = std(rad2deg(z_angle_pp));
z_pp_max = max(rad2deg(z_angle_pp));
z_pp_min = min(rad2deg(z_angle_pp));
omega_pp_mean = mean(omega_pp,2);
omega_pp_sigma = std(omega_pp,0,2);
omega_pp_max = max(omega_pp,[],2);
omega_pp_min = min(omega_pp,[],2);

%Plots
range = max(z_all_max,z_pp_max);
col = 20000;

figure(1)
subplot(2,1,1)
HIST(rad2deg(z_angle_all),col);
ocur1 = HIST(rad2deg(z_angle_all),col);
title('Histogram for NSO Robust Control')
ylabel('Occurrences [.]')
axis([0 0.006 0 max(ocur1)*1.1])
subplot(2,1,2)
HIST(rad2deg(z_angle_pp),col/30)
ocur2 = HIST(rad2deg(z_angle_pp),col/30);
title('Histogram for NSO Simple Control')
ylabel('Occurrences [.]')
xlabel('Error angle [\circ]');
axis([0 0.201 0 max(ocur2)*1.1])
if print_plot
        print(figure(1),'-depsc2', 'main_histogram.eps')
end

% for idx=1:3
% omega_mean_all(idx) = mean(rad2deg(omega_all(idx:end)));
% omega_std_all(idx) = std(rad2deg(omega_all(idx:end)));
% omega_mean_pp(idx) = mean(rad2deg(omega_pp(idx:end)));
% omega_std_pp(idx) = std(rad2deg(omega_pp(idx:end)));
% end

clc
for dddd=1:15
   disp('***********Done!!!!!!***********') 
end
disp('********************** Z-Angle Error **********************');
disp(strcat('NSO structure - Mean:',[32],num2str(z_all_mean),' Std. dev:',[32],num2str(z_all_sigma)))
disp(strcat('Simpel structure - Mean:',[32],num2str(z_pp_mean),' Std. dev:',[32],num2str(z_pp_sigma)))
disp('************************** Omega **************************');
disp(strcat('NSO structure - Mean x:',[32],num2str(omega_all_mean(1)),' Std. dev x:',[32],num2str(omega_all_sigma(1))))
disp(strcat('.....||...... - Mean y:',[32],num2str(omega_all_mean(2)),' Std. dev y:',[32],num2str(omega_all_sigma(2))))
disp(strcat('.....||...... - Mean z:',[32],num2str(omega_all_mean(3)),' Std. dev y:',[32],num2str(omega_all_sigma(3))))
disp(' ')
disp(strcat('Simpel structure - Mean x:',[32],num2str(omega_pp_mean(1)),' Std. dev x:',[32],num2str(omega_pp_sigma(1))))
disp(strcat('......||........ - Mean y:',[32],num2str(omega_pp_mean(2)),' Std. dev y:',[32],num2str(omega_pp_sigma(2))))
disp(strcat('......||........ - Mean z:',[32],num2str(omega_pp_mean(3)),' Std. dev y:',[32],num2str(omega_pp_sigma(3))))
disp(strcat('omega_all_max:',[32],num2str(omega_all_max)))
disp(strcat('omega_all_min:',[32],num2str(omega_all_min)))
disp(strcat('omega_pp_max:',[32],num2str(omega_pp_max)))
disp(strcat('omega_pp_min:',[32],num2str(omega_pp_min)))

%%A method to find the value of the first element in a samle 
%%distribution that is over the value that is the 95% confidence
%%interval value
%sorted(find(sorted >= sorted(length(sorted)*0.95),1,'first'))


%Plot a distribution bell
%sigma = 1/4;
%x1 = -1:1/200:-1/200;
%P = 1.2+1/(sigma*sqrt(2*pi))*exp(-(x1.^2)/(2*sigma^2));
