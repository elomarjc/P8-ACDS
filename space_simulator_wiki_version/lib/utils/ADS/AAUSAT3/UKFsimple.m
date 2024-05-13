function output = UKFsimple(input) 
tic_s=tic;
    % Unscented Kalman filter (UKF) for AAUSAT3. 
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
    % The outputs are all states and residuals
    %
    % Written by Kasper F. Jensen and Kasper Vinther (10gr1035)
    % Aalborg University - 20010
 
    % These variables must be saved between function calls
    persistent x_full;
    persistent P;
    persistent Q;
    persistent R;
    persistent n_xfull;
    persistent n_xerr;
    persistent n_z;
    persistent alpha;
    persistent beta;
    persistent kappa;
    persistent Nrung;
    persistent Ndist;
    persistent I;
    persistent q_s_c;
    persistent lambda;
    persistent n_Xerr;
    persistent W;
    persistent Sun_meas_old;
    persistent Mag_meas_old;
    persistent w_meas_old;
    persistent N_mt_old;
    persistent tictoc;
   
    % Load some of the input data
    t           = input(19);
    q0          = input(20:23);
    Ts          = input(24);
    Eclipse     = input(25);
    
    if t<Ts
  
        % Constants
        n_xfull=7;
        n_xerr=6;
        n_z=9;
        tictoc=[];
        %CoM is 0.049147    0.048996    0.046507 
        %Mass of AAUSAT3: 1.040960 
        %I = diag([0.0018195,0.0024199,0.0024354]);
        %q_s_c=[0.020798,0.7086,-0.019942,0.70502]';
        %q_s_c=q_s_c/norm(q_s_c);
        %CoM is 0.049073    0.048909    0.042976 
        %Mass of AAUSAT3: 0.957820
        I = [0.0017464,0.0022092,0.0022388]'; % Guess of inertia
        q_s_c=[-0.020656,0.71468,-0.007319,0.6991]'; 
        q_s_c=q_s_c/norm(q_s_c);
        alpha=sqrt(3);
        beta=-(1-alpha^2);
        kappa=0;%3-n_xfull;
        Nrung = 10;   % number of RK4 steps during sample time
        Ndist = 0;   % disturbance torques are not modelled on board
        
        % Number of sigma points
        n_Xerr = 2*n_xerr+1;
        
        % Calculate lamda according to scaling parameters
        lambda = alpha^2*(n_xerr+kappa)-n_xerr;
        
        % Array of the weights for each sigma point
        W=[lambda 0.5*ones(1,n_Xerr-1) 0]/(n_xerr+lambda);
        
        % Now calculate the zero'th covariance term weight
        W(n_Xerr+1) = W(1) + (1-alpha^2 + beta);
        
        % Initial covariance matrix:
        vq1=0.001;           % 0.032 standard deviation
        vq2=0.001; 
        vq3=0.001;
        vw1=0.001;           % ca 2 deg/s standard deviation
        vw2=0.001;
        vw3=0.001;
        P=diag([vq1 vq2 vq3 vw1 vw2 vw3]);

        % Model noise matrix:
        vmq1=0.000001;        % 0.001 standard deviation
        vmq2=0.000001;
        vmq3=0.000001;
        vmw1=0.00001;         % ca 0.2 deg/s standard deviation
        vmw2=0.00001; 
        vmw3=0.00001; 
        Q=diag([vmq1 vmq2 vmq3 vmw1 vmw2 vmw3]);

        % Measurement noise matrix
        vSunX=0.0034;
        vSunY=0.0034;
        vSunZ=0.0034;
        vMagX=0.0027; 
        vMagY=0.0027;
        vMagZ=0.0027;
        vwX=0.000012;
        vwY=0.000012;
        vwZ=0.000012;
        R=diag([vSunX vSunY vSunZ vMagX vMagY vMagZ vwX vwY vwZ]);

        % Load omega
        w_meas = input(7:9);

        % Store measurement for next time
        Sun_meas_old = input(1:3);
        Mag_meas_old = input(4:6);
        w_meas_old   = input(7:9);
        N_mt_old     = input(16:18);
     
        % Use Wahba to find q0 if two vector measurements are available
        if Eclipse == 0
            if input(1:3) ~= 0
                if input(4:6) ~= 0
                    % Normalised Sun vector measurements
                    a1=1/vSunX;
                    Sun_I = input(10:12)/(sqrt(input(10:12)'*input(10:12)));
                    Sun_S = input(1:3)/(sqrt(input(1:3)'*input(1:3)));
 
                    % Normalised Magnetic field vector measurements
                    a2=1/vMagX;
                    Mag_I = input(13:15)/sqrt(input(13:15)'*input(13:15));
                    Mag_S = input(4:6)/(sqrt(input(4:6)'*input(4:6)));
 
                    % Use SVD method to solve Wahba's problem
                    [svd_out]=SVD_method(a1,Sun_I,Sun_S,a2,Mag_I,Mag_S,1);
                    q0=svd_out(1:4)';
                    if svd_out(5) > 0.1
                        disp('SVD_method failed to obtain quaternion')
                    end
                end 
            end
        end
                
        % Initial full state
        q0=qmul(q0,q_s_c);
        w_meas(1:4,1)=qmul(qmul(qinv(q_s_c),[w_meas;0]),q_s_c);
        x_full=[q0; w_meas(1:3)];
        
        % Output first time function is run
        omega_s(1:4,1)=qmul(qmul(q_s_c,[x_full(5:7);0]),qinv(q_s_c));
        residual=zeros(9,1);
        output=[qmul(x_full(1:4),qinv(q_s_c)); omega_s(1:3,1); residual];
    else
        
        % Load input data and normalise vectors.
        if input(1:3) ~= Sun_meas_old
            if Eclipse == 0
                Sun_meas= input(1:3)/(sqrt(input(1:3)'*input(1:3)));
                Sun_I   = input(10:12)/(sqrt(input(10:12)'*input(10:12)));
                update_sun = 1;
            else
                update_sun = 0;
            end
        else
            update_sun = 0;
        end
        if input(4:6) ~= Mag_meas_old
            Mag_meas    = input(4:6)/(sqrt(input(4:6)'*input(4:6)));
            Mag_I       = input(13:15)/(sqrt(input(13:15)'*input(13:15)));
            update_mag = 1;
        else
            update_mag = 0;
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
        
        % Calculate sigma points
        X_err = scaledSymmetricSigmaPoints(P,lambda,n_xerr);

        % Expand sigma error states to sigma full states
        X_full=zeros(n_xfull,n_Xerr);
        for ii=1:n_Xerr
            X_full(1:4,ii)=expandstate(X_err(1:3,ii),x_full(1:4));
            X_full(5:7,ii)=x_full(5:7)+X_err(4:6,ii);
        end

        % Predict full state
        X_pfull=zeros(n_xfull,n_Xerr);
        x_pfull=zeros(n_xfull,1);
        N_mt_old(1:4,1)  = qmul(qmul(qinv(q_s_c),[N_mt_old;0]),q_s_c);
        for ii=1:n_Xerr
            X_pfull(:,ii)=runge(Ndist, N_mt_old(1:3,1), X_full(:,ii), I, Ts, Nrung);
            x_pfull=x_pfull+(X_pfull(:,ii)*W(ii)); % mean predicted full state
        end
        x_pfull(1:4)=x_pfull(1:4)./norm(x_pfull(1:4)); %quaternion must be normalized

        % Save torque to next time
        N_mt_old     = input(16:18);
        
        % Convert full state to error state
        X_perr=zeros(n_xerr,n_Xerr);
        x_perr=zeros(n_xerr,1); 
        for ii=1:n_Xerr
            X_perr(1:4,ii)=qmul(X_pfull(1:4,ii),qinv(x_pfull(1:4)));
            X_perr(4:6,ii)=X_pfull(5:7,ii)-x_pfull(5:7);
            x_perr=x_perr+(X_perr(:,ii)*W(ii)); % mean predicted error state
        end
  
        % A priori covariance matrix
        P=zeros(n_xerr);
        for ii=1:n_Xerr
                if ii==1
                    jj=n_Xerr+1;
                else
                    jj=ii;
                end
                P=P+W(jj)*(X_perr(:,ii)-x_perr)*(X_perr(:,ii)-x_perr)';
                % P=P+W(jj)*(X_perr(:,ii))*(X_perr(:,ii))'; %x_perr can be optimized away
        end
        P=P+Q;

        % *************UPDATE****************
        
        % Predict measurements in CRF
        Z_p=zeros(n_z,n_Xerr);
        z_p=zeros(n_z,1);
        for ii=1:n_Xerr
            if update_sun == 1
                Z_p(1:4,ii)=qmul(qmul(qinv(X_pfull(1:4,ii)),[Sun_I;0]),X_pfull(1:4,ii));
            end
            if update_mag == 1
                Z_p(4:7,ii)=qmul(qmul(qinv(X_pfull(1:4,ii)),[Mag_I;0]),X_pfull(1:4,ii));
            end
            if update_w == 1
                Z_p(7:9,ii)=X_pfull(5:7,ii);
            end
            z_p=z_p+(Z_p(:,ii)*W(ii)); % mean predicted measurement
        end
        
        % Rotate measurements from sensors to CRF if new data is available
        % and calculate residual
        residual=zeros(n_z,1); % if there is no update, then residual=0
        if update_sun == 1
            z_meas(1:4,1)=qmul(qmul(qinv(q_s_c),[Sun_meas;0]),q_s_c);
            residual(1:3)=(z_meas(1:3)-z_p(1:3));
        end
        if update_mag == 1
            z_meas(4:7,1)=qmul(qmul(qinv(q_s_c),[Mag_meas;0]),q_s_c);
            residual(4:6)=(z_meas(4:6)-z_p(4:6));
        end
        if update_w == 1
            z_meas(7:10,1)=qmul(qmul(qinv(q_s_c),[w_meas;0]),q_s_c);
            residual(7:9)=(z_meas(7:9)-z_p(7:9));
        end

        % Calculate covariance matrixes based on predicted state and measurements
        P_zz=zeros(n_z,n_z);
        P_xz=zeros(n_xerr,n_z);
        for ii=1:n_Xerr
            if ii==1
                    jj=n_Xerr+1;
            else
                    jj=ii;
            end
            P_zz=P_zz+(W(jj)*(Z_p(:,ii)-z_p)*(Z_p(:,ii)-z_p)');
            P_xz=P_xz+(W(jj)*(X_perr(:,ii)-x_perr)*(Z_p(:,ii)-z_p)');
        end
        P_zz=P_zz+R;

        % Compute Kalman gain
        K=P_xz*P_zz^-1;

        % Update estimated error state
        x_err=K*(residual);

        % Expand error state to full state
        x_full(1:4)=expandstate(x_err(1:3),x_pfull(1:4));
        x_full=[x_full(1:4); x_err(4:6)+x_pfull(5:7)];
                
        % Update covariance matrix
        P=P-(K*P_zz*K');
        eigP=eig(P);
        for j=1:n_xerr
            if eigP(j) <= 0
                error('Error - P is not positive definite after update step')
            end 
        end

        % Output result as a vector
        omega_s(1:4,1)=qmul(qmul(q_s_c,[x_full(5:7);0]),qinv(q_s_c));
        output=[qmul(x_full(1:4),qinv(q_s_c)); omega_s(1:3,1); residual];
  
    end
tictoc=[tictoc,toc(tic_s)];
assignin('base','tictoc_ukf',tictoc);
end
%***************************************************************************

function X = scaledSymmetricSigmaPoints(P,lambda,n_xerr)
    % Find sigma points without mean!
    
    % Calculate matrix square root of weighted covariance matrix
    Psqrtm=(chol((n_xerr+lambda)*P))';
    
    % Array of the sigma points
    X=[zeros(n_xerr,1) -Psqrtm Psqrtm];
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

    temp=(x_err'*x_err);
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

    for kn=1:Nrung
        k1 = f(y       , I, Next);
        k2 = f(y+0.5*stepsize*k1, I, Next);
        k3 = f(y+0.5*stepsize*k2, I, Next);
        k4 = f(y+stepsize*k3    , I, Next);
        y = y+(1/6)*stepsize*(k1+2*k2+2*k3+k4);
        y(1:4)=y(1:4)./norm(y(1:4));
    end
    result = y;
end

%**************************************************************************
function rungres = f(x, I, Next)
    q = x(1:4);
    w = x(5:7);
    I = diag(I);

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
