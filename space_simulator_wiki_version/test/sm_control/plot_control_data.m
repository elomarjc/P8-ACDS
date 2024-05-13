clc             % Clear command window
close all       % Close all figures

plot_time = [280 310];  % Time range for plotting
omega_lim = 1600;       % Limit for angular velocity
u_lim = 3.3;            % Limit for voltage
q_e_lim = [-90,90];     % Limit for quaternion
figure_width = '\figurewidth';  % Figure width for LaTeX
figure_height = '5cm';  % Figure height for LaTeX

reduce_factor = 10;     % Factor for data reduction
LONG_DATA = 0;          % Flag for long data

if LONG_DATA
    for i = 1:4
        q_i_s(:,i) = -1.*decimate(sim_q_s.Data(:,i),reduce_factor,2);   % Decimate satellite attitude quaternion
        q_i_ref(:,i) = decimate(sim_q_ref.Data(:,i),reduce_factor,2);     % Decimate reference quaternion
    end
    q_i_s(:,5) = decimate(sim_q_s.Time,reduce_factor,2);   % Decimate time
    q_i_ref(:,5) = decimate(sim_q_ref.Time,reduce_factor,2); % Decimate time
    
    % Plot attitude quaternions
    figure;
    plot(q_i_s(:,5),q_i_s(:,1),q_i_s(:,5),q_i_s(:,2),q_i_s(:,5),q_i_s(:,3),...
        q_i_s(:,5),q_i_s(:,4));
    hold on
    plot(q_i_ref(:,5),q_i_ref(:,1),'k--',q_i_ref(:,5),q_i_ref(:,2),'k--',...
        q_i_ref(:,5),q_i_ref(:,3),'k--',q_i_ref(:,5),q_i_ref(:,4),'k--');
    xlim([200 800]);
    title('')
    ylabel('Quaternion [-]');
    xlabel('Time [s]');
    legend('q_{1}','q_{2}','q_{3}','q_{4}');
    grid on
    grid minor
    matlab2tikz('fig/attitude_long.tex',...
        'height',figure_height,...
        'width',figure_width)
else
    q_i_s = sim_q_s;
    q_i_ref = sim_q_ref;
    
    % Plot attitude quaternions
    figure;
    plot(-1*q_i_s);
    hold on
    p_ref = plot(q_i_ref,'k--');
    xlim(plot_time);
    title('')
    ylabel('Quaternion [-]');
    xlabel('Time [s]');
    legend('q_{1}','q_{2}','q_{3}','q_{4}');
    grid on
    grid minor
    matlab2tikz('fig/attitude.tex','height',figure_height,'width',figure_width)
end

euler = sim_euler_e;    % Euler angle error
u_rw = sim_u_rw;        % Reaction wheel voltage
omega_rw = sim_omega_rw;    % Reaction wheel angular velocity

% Plot Euler angle error
figure;
p_e = plot(euler);
hold on
plot(plot_time,[5 5],'k--',plot_time,[-5 -5],'k--')
xlim(plot_time);
ylim(q_e_lim);
title('')
ylabel('Angle Error [deg]');
xlabel('Time [s]');
legend('Roll','Pitch','Yaw');
grid on
grid minor
matlab2tikz('fig/euler_error.tex','height',figure_height,'width',figure_width)

% Plot reaction wheel voltage
figure;
p_u = plot(u_rw);
xlim(plot_time);
ylim([0,u_lim]);
title('')
ylabel('Voltage [V]');
xlabel('Time [s]');
legend('RW1','RW2','RW3','RW4');
grid on
grid minor
matlab2tikz('fig/control.tex','height',figure_height,'width',figure_width)

% Plot reaction wheel angular velocity
figure;
p_omega = plot(omega_rw);
xlim(plot_time);
ylim([0,omega_lim]);
title('')
ylabel('Angular velocity [rad/s]');
xlabel('Time [s]');
legend('RW1','RW2','RW3','RW4');
grid on
grid minor
matlab2tikz('fig/motor_vel.tex','height',figure_height,'width',figure_width)
