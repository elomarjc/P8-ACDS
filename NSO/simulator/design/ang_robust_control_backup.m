%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the design file for the robust angular velocity controller.
% It requires YALMIP with SeDuMi and calculates the state feedback
% controller complying with the strongly robust H_inf performance criterion.
% Furthermore, a pole placement restraint is introduced.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
analyse = 1;  %% Turn this ON for plots.
one_plot = true;
plot_pole_zero = true;
system_ang = 2; % 1: 0.5(w - w_r) 2: no w3
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

% System description
A_sys = [inv(I)*Sh];
B_sys = -inv(I);
C_sys = eye(3);

sys = ss(A_sys,B_sys,C_sys,0);

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
if system_ang == 1 % system with weight on (0.5w - 0.5w_r)
    for k_ang=1:N
        I_delta = diag(inv(I)*corners(k_ang,:)');
        A{k_ang}  =  I_delta*Sh;
        B1{k_ang} =  1/gamma*[I_delta zeros(3)];
        B2{k_ang} = -I_delta;
    end
    D12 = zeros(3);
    C1  = [0.5*eye(3)];
    D11 = [zeros(3) -0.5*eye(3)];
elseif system_ang == 2 % system with no w3 in z
    for k_ang=1:N
        I_delta = diag(inv(I)*corners(k_ang,:)');
        A{k_ang}  =  I_delta*Sh;
        B1{k_ang} =  1/gamma*[I_delta];
        B2{k_ang} = -I_delta;
    end
    D12 = zeros(2,3);
    C1  = [eye(2) zeros(2,1)];
    D11 = zeros(2,3);
elseif system_ang == 3 % system for P8-ACDS
    for k_ang=1:N
        I_delta = diag(inv(I)*corners(k_ang,:)');
        A{k_ang}  =  I_delta*Sh;
        B1{k_ang} =  1/gamma*[I_delta];
        B2{k_ang} = -I_delta;
    end
    D12 = zeros(2,3);
    C1  = [eye(2) zeros(2,1)];
    D11 = zeros(2,3);
end

% dimensions used in the LMI construction
n = size(A{1},1);
q = size(B1{1},2);
m = size(B2{1},2);
n1 = size(D12,1);

yalmip('clear');

% Define variables
W = sdpvar(m, n);
Y = sdpvar(n, n, 'symmetric');
theta = 20; % Change the value of theta as per your requirement
epsilon = 1e-6; % Small positive value

% Define R1
R1 = eye(q) - D11'*D11;

% Define the additional constraint matrix
C = [sin(theta) * ((Y * A{k_ang} + W * B2{k_ang} + (A{k_ang}') * Y + (B2{k_ang}') * W')), ...
     cos(theta) * ((Y * A{k_ang} + W * B2{k_ang} - (A{k_ang}') * Y + (B2{k_ang}') * W'));
     cos(theta) * (((A{k_ang}') * Y + (B2{k_ang}') * W' + Y * A{k_ang} + W * B2{k_ang})), ...
     sin(theta) * ((Y * A{k_ang} + W * B2{k_ang} + (A{k_ang}') * Y + (B2{k_ang}') * W'))];

% Define the constraint using a small margin to ensure numerical satisfaction
Constraints = [Y >= epsilon * eye(n), C <= -epsilon];

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
            plot(pole(feedback(F_ang*ss(A{k_ang}, B2{k_ang}, eye(3), 0), eye(3), 1)), 'xb')
        end
        grid on
        if one_plot
            hold on
            plot(pole_lim*ones(2,1), 0.2*[1 -1], 'r')
            plot(pole(feedback(F_ang*ss(A_sys, B_sys, eye(3), 0), eye(3), 1)), '^k')
        end
    end % plot poles
end % end analyse==1

if print_plot && one_plot
    print -depsc2 'ang_control_design_poles.eps'
end

% Used to calculate the pole placement feedback controller
P_ang = [-1.25 -1.26 -1.27]'; 
K_ang = -place(A_sys, B_sys, P_ang)   %% Show the K_ang

if clear_variables
    clear K_ang P_ang F_ang W Y A A_sys B1 B2 B_sys C1 C_sys D11 D12 F11 F12 F13 FP I I_delta ...
          N Phi R1 Sh W_feas Y_feas analyse clsys corners gamma h hn k_ang lmiset ...
          m m_do m_up n n1 one_plot options plot_pole_zero pole_lim print_plot ...
          q solution sys system_ang clear_variables run_silent
end
