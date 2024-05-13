%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the design file for the robust angular velocity controller.
% It requires YALMIP with SeDuMi and calculates the state feedback
% controller complying with the strongly robust H_inf performance criterion.
% Furthermore, a pole placement restraint is introduced.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
analyse = 1;  %% Turn this ON for plots.
one_plot = true;
plot_pole_zero = true;
system_ang = 1; % 1: P8-adcs
print_plot = false;
clear_variables = true;
run_silent = true;

% operating point
operating_point_constants;

% limits on inertia matrix
m_up = 1.05;
m_do = 0.95;

% pole limit
pole_lim = -0.2*2*pi;

% Dummy variables
w_0 = [0, 0, 0.1]';
q_0 = [0, 0, 0]';

% Inertia matrix
I_p = [0.0088, 0, 0; 
       0, 0.0088, 0; 
       0, 0, 0.0044];

% System description
A_sys = [[inv(I_p)*(skew3(I_p*w_0)-skew3(w_0)*I_p), zeros(3);
         -1/2*(skew3(q_0)+eye(3))                         , 1/2*skew3(w_0)],zeros(6,3);
         zeros(3),eye(3),zeros(3)];
B_sys =[-inv(I_p);zeros(6,3)];
C_sys = [eye(9)];

sys = ss(A_sys,B_sys,C_sys,0);

% A_sys = [inv(I_p)*(skew3(I_p*w_0)-skew3(w_0)*I_p), zeros(3);
%          -1/2*skew3(q_0)                         , 1/2*skew3(w_0)];
%B_sys = [-inv(I_p);zeros(3)];
%C_sys = [eye(6)];

% dimensions used in the LMI construction
n = 9;
q = 3;
m = 3;
n1 = 3;

yalmip('clear');

% Define variables
W = sdpvar(m, n);
Y = sdpvar(n, n, 'symmetric');
theta = 20; % Change the value of theta as per your requirement
epsilon = 1e-6; % Small positive value

% vertices of the convex uncertainty region
corners = [m_do m_do m_do;
          m_up m_do m_do;
          m_up m_up m_do;
          m_do m_up m_do;
          m_do m_do m_up;
          m_up m_do m_up;
          m_up m_up m_up;
          m_do m_up m_up];

% number of corners
N = size(corners,1);

gamma = 20;
if system_ang == 1 % system for P8-ACDS
    for k_ang=1:N
        I_delta = diag(inv(I)*corners(k_ang,:)');
        A{k_ang}  =  [[inv(I_delta)*(skew3(I_delta*w_0)-skew3(w_0)*I_delta), zeros(3);
         -1/2*(skew3(q_0)+eye(3))                         , 1/2*skew3(w_0)],zeros(6,3);
         zeros(3),eye(3),zeros(3)];
        B1{k_ang} = ones(9,3);
        B2{k_ang} = [-inv(I_delta);zeros(6,3)];
    end
    D12 = zeros(3,3);
    C1  = [eye(3) zeros(3,6)];
    D11 = -eye(3,3);
end

% Define R1
R1 = eye(q) - D11'*D11;

% Define the additional constraint matrix
% C = [sin(theta) * (Y * A{k_ang} + W * B2{k_ang} + A{k_ang}' * Y + B2{k_ang}' * W'), ...
%      cos(theta) * (Y * A{k_ang} + W * B2{k_ang} - A{k_ang}' * Y + B2{k_ang}' * W');
%      cos(theta) * (A{k_ang}' * Y + B2{k_ang}' * W' + Y * A{k_ang} + W * B2{k_ang}), ...
%      sin(theta) * (Y * A{k_ang} + W * B2{k_ang} + A{k_ang}' * Y + B2{k_ang}' * W')];

% Define the constraint using a small margin to ensure numerical satisfaction
%Constraints = [Y >= epsilon * eye(n), C <= -epsilon];
Constraints = [Y >= epsilon * eye(n)];

for k_ang=1:N
    F11{k_ang} = Y*A{k_ang}' + A{k_ang}*Y + W'*B2{k_ang}' + B2{k_ang}*W;
    F12{k_ang} = B1{k_ang} + Y*C1'*D11 + W'*D12'*D11;
    F13{k_ang} = Y*C1' + W'*D12';

    Phi{k_ang} = [F11{k_ang}, F12{k_ang}, F13{k_ang};
                  F12{k_ang}', -R1, zeros(q, n1);
                  F13{k_ang}', zeros(n1, q), -eye(n1)];

    % Pole placement constraint
    FP{k_ang} = Y*2*pole_lim*eye(n) - Y*A{k_ang}' - W'*B2{k_ang}' - A{k_ang}*Y - B2{k_ang}*W;

    % Add LMIs to the constraint list
    Constraints = [Constraints, FP{k_ang} <= 0, Phi{k_ang} <= 0];
end

% Set YALMIP options
options = sdpsettings('solver', 'sedumi', 'verbose', ~run_silent);
%options = sdpsettings('solver', 'mosek', 'verbose', ~run_silent, 'sedumi.eps', 1e-6);

% Solve the problem
solution = optimize(Constraints, trace(Y), options);

% Retrieve the values from the YALMIP variable
Y_feas = value(Y);
W_feas = value(W);

% Check the solution and display messages if necessary
if ~run_silent
    check(Constraints);
end

if solution.problem
    disp('******** solver could not find a solution ********')
end

% Calculate the controller
F_ang = W_feas / Y_feas;

% Define the variables used in the angular controller block in Simulink
A_ang_c = zeros(3, 3);
B_ang_c = zeros(3, 3);
C_ang_c = zeros(3, 3); 
D_ang_c = F_ang;

% fprintf('Size of A_sys: %dx%d\n', size(A_sys));
% fprintf('Size of B_sys: %dx%d\n', size(B_sys));
% fprintf('Size of F_ang: %dx%d\n', size(F_ang));
% fprintf('Size of C1: %dx%d\n', size(C1));
% fprintf('Size of D12: %dx%d\n', size(D12));
% fprintf('Size of D11: %dx%d\n', size(D11));

% Used to analyze the designed controller
if analyse
    for k_ang=1:N
        clsys = pck(A{k_ang} + B2{k_ang}*F_ang, B1{k_ang}, C1 + D12*F_ang, D11);
        hn = hinfnorm(clsys);
        %disp(hn(2)); % why display Hn2?
    end
    clsys = pck(A_sys + B_sys*F_ang, B_sys/gamma, C1 + D12*F_ang, D11);
    hn = hinfnorm(clsys);
    % disp(hn(2));  % why display Hn2?
    
    if plot_pole_zero
        if one_plot
            figure(1)
            clf
            hold on
        end
        for k_ang=1:N
            if ~one_plot
                figure(k_ang)
                clf
                hold on
                plot(pole_lim*ones(2,1), 0.2*[1 -1], 'r')
                grid on
            end
            plot(pole(feedback(F_ang*ss(A{k_ang}, B2{k_ang}, eye(9), 0), eye(3), 1)), 'xb')
        end
        grid on
        if one_plot
            hold on
            plot(pole_lim*ones(2,1), 0.2*[1 -1], 'r')
            plot(pole(feedback(F_ang*ss(A_sys, B_sys, eye(9), 0), eye(3), 1)), '^k')
        end
    end % plot poles
end % end analyse==1

if print_plot && one_plot
    print -depsc2 'ang_control_design_poles.eps'
end

% Used to calculate the pole placement feedback controller
%P_ang = [-1.25, -1.28, -1.31, -1.22, -1.25, -1.28, -1.19, -1.22, -1.25]';
%K_ang = -place(A_sys, B_sys, P_ang)   %% Show the K_ang

if clear_variables
    clear K_ang P_ang F_ang W Y A A_sys B1 B2 B_sys C1 C_sys D11 D12 F11 F12 F13 FP I I_delta ...
          N Phi R1 Sh W_feas Y_feas analyse clsys corners gamma h hn k_ang lmiset ...
          m m_do m_up n n1 one_plot options plot_pole_zero pole_lim print_plot ...
          q solution sys system_ang clear_variables run_silent
end

function [skew]=skew3(u)
    skew = cross(repmat(u,1,3),eye(3));
end
