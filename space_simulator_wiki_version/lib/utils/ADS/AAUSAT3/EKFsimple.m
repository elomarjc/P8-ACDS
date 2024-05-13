function output = EKFsimple(input)
% tic
% Extended Kalman filter (EKF) for AAUSAT3. 
%
% The filter states are:
% x=[q_i_s omega_s]'=[q1 q2 q3 q4 w1 w2 w3]'
%
% The input measurements are a vector to the Sun (SBRF), 
% a magnetic field vector (SBRF) and omega (SBRF).
% The filter also needs a calculated Sun vector (ECI), 
% a calculated magnetic field vector (ECI), eclipse 
% indication and an initial quaternion (q_i_s). 
%
% Written by Kasper F. Jensen and Kasper Vinther (09gr935)
% Aalborg University - 2009

% These variables must be saved between function calls
    persistent x_full;
    persistent P;
    persistent Q;
    persistent R;
    persistent n_xerr;
    persistent n_z;
    persistent Nrung;
    persistent Ndist;
    persistent I;
    persistent q_s_c;
    persistent Sun_meas_old;
    persistent Mag_meas_old;
    persistent w_meas_old;
    persistent N_mt_old;
    
    % Load some of the input data
    t           = input(19);
    q0          = input(20:23);
    Ts          = input(24);
    Eclipse     = input(25);
    
    if t < Ts

        % Constants
        n_xerr=6;
        n_z=9;
        %CoM is 0.049073    0.048909    0.042976 
        %Mass of AAUSAT3: 0.957820
        I = [0.0017464,0.0022092,0.0022388]'; % Guess of inertia
        q_s_c=[-0.020656,0.71468,-0.007319,0.6991]'; 
        q_s_c=q_s_c/norm(q_s_c);
        %I = diag([0.003,0.003,0.003]);
        %q_s_c=[0 0 0 1]';
        
        Nrung = 20;   % number of steps during sample time
        Ndist = 0;   % disturbances are not modelled on board
        
        % Initial covariance matrix:
        vq1=0.01;
        vq2=0.01;
        vq3=0.01;
        vw1=0.001;
        vw2=0.001;
        vw3=0.001;
        P=diag([vq1 vq2 vq3 vw1 vw2 vw3]);

        % Model noise matrix:
        vmq1=0.0001;
        vmq2=0.0001;
        vmq3=0.0001;
        vmw1=0.00001;
        vmw2=0.00001;
        vmw3=0.00001;
        Q=diag([vmq1 vmq2 vmq3 vmw1 vmw2 vmw3]);

        % Measurement noise matrix
        vSunX=0.031;
        vSunY=0.031;
        vSunZ=0.031;
        vMagX=0.0012; 
        vMagY=0.0012;
        vMagZ=0.0012;
        vwX=0.000012;
        vwY=0.000012;
        vwZ=0.000012;
        R=diag([vSunX vSunY vSunZ vMagX vMagY vMagZ vwX vwY vwZ])/10;
        
        % Load omega
        w_meas = input(7:9);

        % Store measurement for next time
        Sun_meas_old = input(1:3);
        Mag_meas_old = input(4:6);
        w_meas_old   = input(7:9);
        N_mt_old     = input(16:18);
        
        % Initial full state
        q0=qmul(q0,q_s_c);
        w_meas(1:4,1)=qmul(qmul(qinv(q_s_c),[w_meas;0]),q_s_c);
        x_full=[q0; w_meas(1:3)];
        
        % Output first time function is run
        omega_s(1:4,1)=qmul(qmul(q_s_c,[x_full(5:7);0]),qinv(q_s_c));
        output=[qmul(x_full(1:4),qinv(q_s_c)); omega_s(1:3,1)];
    else

        % Load input data and normalise vectors.
        if input(1:3) ~= Sun_meas_old
            if Eclipse == 0
                Sun_meas= input(1:3)/(sqrt(input(1:3)'*input(1:3)));
                Sun_I   = input(10:12)/(sqrt(input(10:12)'*input(10:12)));
                update_sun = 1;
            else
                update_sun = 0;
                %   Sun_I = [0 0 0]';
            end
        else
            update_sun = 0;
            % Sun_I = [0 0 0]';
        end
        if input(4:6) ~= Mag_meas_old
            Mag_meas    = input(4:6)/(sqrt(input(4:6)'*input(4:6)));
            Mag_I       = input(13:15)/(sqrt(input(13:15)'*input(13:15)));
            update_mag = 1;
        else
            update_mag = 0;
            %  Mag_I = [0 0 0]';
        end
        if input(7:9) ~= w_meas_old
            w_meas      = input(7:9);
            update_w = 1;
        else
            update_w = 0;
        end

        % Store measurement for next time
        Sun_meas_old = input(1:3);
        Mag_meas_old = input(4:6);
        w_meas_old   = input(7:9);
        
        % ************PREDICT****************
        
        % Predict full state
        N_mt_old(1:4,1)  = qmul(qmul(qinv(q_s_c),[N_mt_old;0]),q_s_c);
        x_pfull = runge(Ndist, N_mt_old(1:3,1), x_full, I, Ts, Nrung);
        N_mt_old     = input(16:18);

        % Calculate Jacobian matrix
        Phi = jacobianF(x_full,I,Ts); 
        %eigenvalus_PHIsys=eig(Phi)
        
        % Propagate the covariance matrix:
        P_before=P;
        P   = Phi*P*Phi' + Q;
        
        % *************UPDATE****************
        
        % Predict measurements in CRF and rotate measurements from sensors 
        % to CRF if new data is available
        z_p=zeros(n_z,1);
        if update_sun == 1
            z_p(1:4,1)=qmul(qmul(qinv(x_pfull(1:4)),[Sun_I;0]),x_pfull(1:4));
            z_meas(1:4,1)=qmul(qmul(qinv(q_s_c),[Sun_meas;0]),q_s_c);
        else
            z_meas(1:3,1)=z_p(1:3)
        end
        if update_mag == 1
            z_p(4:7,1)=qmul(qmul(qinv(x_pfull(1:4)),[Mag_I;0]),x_pfull(1:4));
            z_meas(4:7,1)=qmul(qmul(qinv(q_s_c),[Mag_meas;0]),q_s_c);
        else
            z_meas(4:6,1)=z_p(4:6);
        end
        if update_w == 1
            z_p(7:9,1)=x_pfull(5:7);
            z_meas(7:10,1)=qmul(qmul(qinv(q_s_c),[w_meas;0]),q_s_c);
            z_meas=z_meas(1:9);
        else
            z_meas(7:9,1)=z_p(7:9);
        end

        % Calculate Jacobian measurement matrix
        %z_sun_bar(1:4,1)=qmul(qmul(qinv(x_full(1:4)),[Sun_I;0]),x_full(1:4));
        %z_mag_bar(1:4,1)=qmul(qmul(qinv(x_full(1:4)),[Mag_I;0]),x_full(1:4));
        
        H = [2*skew_matrix(z_p(1:3)) zeros(3);
             2*skew_matrix(z_p(4:6)) zeros(3);
             zeros(3)                   eye(3) ];
        
        % Compute Kalman gain
        K = P*H'*((H*P*H' + R)^(-1));

        % Calculate error state
        x_err=K*(z_meas-z_p)
        %eigenvalues=eig(Phi-K*H)
        
        % Expand error state to full state
        x_full(1:4)=expandstate(x_err(1:3),x_pfull(1:4));
        x_full=[x_full(1:4); x_err(4:6)+x_pfull(5:7)];

        % Update covariance matrix
        P = (eye(n_xerr) - K*H)*P;
        %P_diff=P-P_before
        
        % Output result as a vector
        omega_s(1:4,1)=qmul(qmul(q_s_c,[x_full(5:7);0]),qinv(q_s_c));
        output=[qmul(x_full(1:4),qinv(q_s_c)); omega_s(1:3,1)];
    end
    %  toc
end

%***************************************************************************
function PHIsys = jacobianF(x,I,Ts)
% Returns the Jacobian matrix of the model

    w       = x(5:7);
    
    S_w     = skew_matrix(w);
    S_Iw    = skew_matrix(I*w);
    
    Fsys    = [ -S_w        .5*eye(3)         ;
                zeros(3)    I^(-1)*(S_Iw-S_w*I) ];
    %eigenvalus_con=eig(Fsys)
    %Gsys    = [zeros(3) ; 
    %            I^(-1)  ]; 
    PHIsys = eye(6)+Ts*Fsys;
    %eigenvalus_PHIsys2=eig(PHIsys2)
    %ss_sys  = ss(Fsys,Gsys,eye(6),0);
    %dss_sys = c2d(ss_sys,Ts);                
    %PHIsys  = dss_sys.A
end

%**************************************************************************
function output = skew_matrix(x)
% Returns a skew symmetric matrix of the input vector.

    if length(x) ~= 3
        error('skew_matix: dimension error, - input is not a 1x3 vector')
    end
    
    output =  [0    -x(3)  x(2);
               x(3)  0    -x(1);
               -x(2)  x(1)  0  ];
end

%**************************************************************************
function x_exp = expandstate(x_err,x_full)
% Expands the quaternion state from 3 elements to 4
    x_err
    temp=(x_err'*x_err)
    if temp<=1
        x_exp=[x_err; sqrt(1-temp)];
    else
        x_exp=[x_err/sqrt(1+temp); 1/sqrt(1+temp)]; % [adcs_04lars,76 & Humphreys,2002]
        disp('Filter in trouble - x_err unusually large')
    end
    x_exp=qmul(x_exp,x_full);
end

%**************************************************************************
function result = runge(Ndist, Nmt, x0, I, step, Nrung)
% Step is the time to jump forward
% Nrung is the number of runs to do it in

    Next = Ndist + Nmt;

    y=x0;
    stepsize = step/Nrung;

    for kn=1:Nrung,
        k1 = f(y       , I, Next);
        k2 = f(y+0.5*stepsize*k1, I, Next);
        k3 = f(y+0.5*stepsize*k2, I, Next);
        k4 = f(y+stepsize*k3    , I, Next);
        y = y+(1/6)*stepsize*(k1+2*k2+2*k3+k4);
        y(1:4)=y(1:4)./sqrt(y(1:4)'*y(1:4));
    end

    result = y;
end

%**************************************************************************
function rungres = f(x, I, Next)
    q = x(1:4);
    w = x(5:7);

    q_dot = 0.5.*[-skew_matrix(w) w; -w' 0]*q;
    w_dot = I^(-1)*(-skew_matrix(w)*(I*w) + Next);

    rungres = [q_dot; w_dot];
end

%**************************************************************************
function y = qinv(q)
    y = [-q(1:3); q(4)];
end

%**************************************************************************
function qres = qmul(q1, q2)
    qres1 =  q1(1)*q2(4) + q1(2)*q2(3) - q1(3)*q2(2) + q1(4)*q2(1);
    qres2 = -q1(1)*q2(3) + q1(2)*q2(4) + q1(3)*q2(1) + q1(4)*q2(2);
    qres3 =  q1(1)*q2(2) - q1(2)*q2(1) + q1(3)*q2(4) + q1(4)*q2(3);
    qres4 =  q1(4)*q2(4) - q1(1)*q2(1) - q1(2)*q2(2) - q1(3)*q2(3);

    qres = [qres1 qres2 qres3 qres4]';
end
