analyse = 1;
% operating point
h     = [2.23 2.23 2.23]*1e-3;       % momentum wheel bias
Sh    = [0 -h(3) h(2);h(3) 0 -h(1);-h(2) h(1) 0];
I     = diag([42.3 42.3 28.4]*1e-3); % inertia matrix for deployed situation

m_up = 1.01;
m_do = 0.99;


% System description
A_sys = [inv(I)*Sh];
B_sys = -inv(I);
C_sys = eye(3);

sys = ss(A_sys,B_sys,C_sys,0);

corners = [m_do m_do m_do;
          m_up m_do m_do;
          m_up m_up m_do;
          m_do m_up m_do;
          m_do m_do m_up;
          m_up m_do m_up;
          m_up m_up m_up;
          m_do m_up m_up];

N = size(corners,1);

epsilon = 1;
gamma = 1/5;
k = 1;
for l=1:N
    I_delta = diag(inv(I)*corners(l,:)');
    A{l}  =  I_delta*Sh;
    B1{l} =  [(1/gamma)*I_delta zeros(3)];
    B2{l} = -I_delta*k;
end

D12=  [zeros(2,3);epsilon*eye(3)];
C1  = [0.5*[eye(2) zeros(2,1)];zeros(3)];%[eye(3);zeros(3)];%
D11 = [zeros(2,3) -0.5*[eye(2) zeros(2,1)];zeros(3) zeros(3)];%zeros(6,3);%

n = size(A{1},1);
q = size(B1{l},2);
m = size(B2{1},2);
n1= size(D12,1);

yalmip('clear');

W = sdpvar(m,n);
Y = sdpvar(n,n,'symmetric');

R1 = eye(q)-D11'*D11;

% Define the constraint using a semidefinite constraint
Constraints = [Y >= 0]; 

for l=1:N
    F11{l} = Y*A{l}'+A{l}*Y+W'*B2{l}' + B2{l}*W;
    F12{l} = B1{l} +Y*C1'*D11 + W'*D12'*D11;
    F13{l} = Y*C1'+W'*D12';
  
    Phi{l} = [F11{l}   F12{l}      F13{l};
            F12{l}'  -R1       zeros(q,n1);
            F13{l}' zeros(n1,q)   -eye(n1)];
  
    Constraints = [Constraints, Phi{l} <= 0];
end

options = sdpsettings('solver','sedumi'); %sdpt3
solution = optimize(Constraints, trace(Y), options);
if solution.problem
    disp('******** YALMIP could not find a solution ********')
end

Y_feas = value(Y);
W_feas = value(W);

% Check the feasibility of constraints
check(Constraints)

F = W_feas/Y_feas;

A_ang_c = zeros(3,3);
B_ang_c = zeros(3,3);
C_ang_c = zeros(3,3); 
D_ang_c = F;

if analyse
    NoP = 50;
    w = logspace(-1,3,NoP);
    for l=1:N
        clsys = pck(A{l}+B2{l}*F,B1{l},C1+D12*F,D11);
        hn = hinfnorm(clsys);
        disp(hn(2));
        sys_cl{l} = ss(A{l}+B2{l}*F,B1{l},C1+D12*F,D11);
        mag{l} = bode(sys_cl{l},w);
        smag{l} = squeeze(mag{l});
        lsv{l} = sqrt(sum(abs(smag{l}).^2, 1)); %lsv{l} = sqrt(diag(smag{l}'*smag{l}));
        
        % Debugging statements
        disp(['Length of lsv{l}: ', num2str(length(lsv{l}))]);
        disp(['Length of w: ', num2str(length(w))]);
        disp(['l: ', num2str(l)]);
        disp(['N: ', num2str(N)]);
        
        figure;
        loglog(w, reshape(lsv{l}, [], 1), w, ones(NoP,1));
        title(['Magnitude plot for corner ', num2str(l)]);
        xlabel('Frequency (rad/s)');
        ylabel('Magnitude');
        grid on;
        legend('Magnitude', 'Unity Magnitude');
    end
end

% Desired poles
P=[-1.25 -1.26 -1.27]'; 

K=-place(A_sys,B_sys,P);
