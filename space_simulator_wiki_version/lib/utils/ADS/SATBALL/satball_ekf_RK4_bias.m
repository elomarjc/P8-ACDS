function [x_hat,phi,H, P_apo_t] = satball_ekf_RK4_bias(z_meas, u, x0, q_s_c, q_s_tetra,t)
    %> @brief Multiplicative Extended Kalman Filter for use in the attitude determination 
    %> @brief of the Aalborg University Studentspace Attitude Control Teststand.
    %>
    %> Use:
    %> [x_hat,phi,H] = ekf(z_meas, u, x0, dt, q_s_c, q_s_tetra)
    %> 
    %> Subfunctions:
    %>   ekf_prediction()
    %>   ekf_update()
    %>   ekf_jacobian()
    %>   reaction_wheels()
    %>
    %> @param z_meas: measurement vector containing (acc, gyro, mag, rw]
    %> @param u: control inputs to the reaction wheels
    %> @param x0: initial statevector of the system
    %> @param dt: sampletime
    %> @param q_s_c: quaternion describing the rotation from body to control frame
    %> @param q_s_tetra: quaternion describing the rotation from body to tetrahedron
    %>        
    %> @retval x_hat: aposteriori state vector estimate [q, w, rw]
    %> @retval phi: system matrix
    %> @retval H: measurement matrix

    assert(size(z_meas,1) == 13, 'satball_ekf: z_meas is wrong dimension');
    assert(size(u,1) == 4, 'satball_ekf: u is wrong dimension');
    assert(size(x0,1) == 14, 'satball_ekf: x0 is wrong dimension');
    assert(size(q_s_c,1) == 4, 'satball_ekf: q_s_c is wrong dimension');
    assert(size(q_s_tetra,1) == 4, 'satball_ekf: q_s_tetra is wrong dimension');

    persistent P_apo x_apo first J Ji Q R tstart t_old;

    if isempty(first)
        % initialisation

        % Process noise covariance matrix
        q = [ones(1,6)*1e-6,...  % Model states [q_bar_s, omega_s]'
             ones(1,3)*1e-10,... % gyro bias covariance
             ones(1,4)*1e-3];    % Model reaction wheels [rw_omega]

        Q = diag(q);
        
        p = [ones(1,6)*1e-1,...  % Model states [q_bar_s, omega_s]'
             ones(1,3)*1e-8,...  % gyro bias covariance
             ones(1,4)*1e-1];    % Model reaction wheels [rw_omega]
        P_apo = diag(p);

        % Measurement noise covariance matrix
        r = [ones(1,3)*1e-1,...    % Accelerometer mesurement covariance
             ones(1,3)*3.4e-6,...  % Gyro covariance
             ones(1,3)*2.23e-3,... % Magnetometer covariance
             ones(1,4)*1e-1];      % Reaction wheel covariance            
        
        R = diag(r);

        x_apo = zeros(14,1);
        
        J = diag([0.6714 0.7044 0.7268]*1e-3);
        Ji = pinv(J);
        t_old = t;
        first = 1;
    end

    % calculate sample time
    Ts = t - t_old;
    t_old = t;

    % Prediction step
    if first == 1
        [x_apri, z_apri, P_apri, H, phi] = ekf_prediction(x0, u, J, Ji, Ts, Q, P_apo, q_s_c, q_s_tetra);
        first = 0;
    else
        [x_apri, z_apri, P_apri, H, phi] = ekf_prediction(x_apo, u, J, Ji, Ts, Q, P_apo, q_s_c, q_s_tetra);
    end
    
    % Update step
    [x_apo_c, P_apo] = ekf_update(x_apri, z_meas, z_apri, H, P_apri, R, q_s_c);
    
    % rotating and returning state vector into satellite frame
    % quaternion
    x_apo(1:4,1) = qmult(x_apo_c(1:4), qinv(q_s_c)); 
    % omega
    x_apo(5:7,1) = qRot(x_apo_c(5:7), qinv(q_s_c));
    % gyro_bias
    x_apo(8:10,1) = qRot(x_apo_c(8:10,1), qinv(q_s_c));
    % Reaction Wheels
    x_apo(11:14,1) = x_apo_c(11:14);
    
    x_hat = x_apo;
    P_apo_t = P_apo;
end


% ======================================================================
%> @brief Running the prediction step of the Extended Kalman Filter
%>
%> @param x_apo: Previousaposteriori state estimate [q, w, rw, ab, gb, mb]
%> @param u: Current actuator input
%> @param J: Inertia matrix of the system
%> @param Ji: Inverse inertial matrix
%> @param dt: Timestep
%> @param Q:
%> @param P_apo:
%> @param q_s_c:
%> @param q_s_tetra:
%>
%> @retval x_apri: Apriori state vector estimate
%> @retval z_apri: Apriori measurement vector estimate
%> @retval P_apri:
%> @retval H:
%> @retval phi:
% ======================================================================
function [x_apri, z_apri, P_apri, H, phi] = ekf_prediction(x_apo, u, J, Ji, dt, Q, P_apo, q_s_c, q_s_tetra)
    % Checking for dimensions of the input arguments
    assert(size(x_apo,1) == 14, 'ekf_prediction: x_apo wrong dimension');
    assert(size(u,1) == 4, 'ekf_prediction: u wrong dimension');
    
    % Gravity
    g = [0 0 1]';
    mag_t = [0 1 0]';
    % mag_t = [0.0027    0.9817   -0.1903]';
    
    rk_rate = 20;

    % Separating the state vector into variables and rotating to control frame
    % Satellite Attitude
    q_apo_b = x_apo(1:4);
    q_apo_c = qmult(q_apo_b, q_s_c);

    % Angular velocity
    w_apo_b = x_apo(5:7);
    w_apo_c = qRot(w_apo_b, q_s_c);
    
    % gyro bias
    gyro_bias_c = qRot(x_apo(8:10), q_s_c); % w_dot_bias = 0 + v  =>  apri = apo
    
    % Reaction Wheel Velocity
    rw_omega_apo = x_apo(11:14);

    % Actuator torque and coriolis
    act_torque = reaction_wheels(q_s_c, q_s_tetra, w_apo_b, rw_omega_apo, u);

    % State Propagation using 4th order Runge-Kutta method
    x_apri_c = rk4(J, Ji, w_apo_c, act_torque, q_apo_c, rw_omega_apo, u, dt, rk_rate);

    q_apri_c = x_apri_c(1:4);
    w_apri_c = x_apri_c(5:7);
    rw_omega_apri = x_apri_c(8:11);

    % Mapping the apriori state estimate to the apriori measurement vector
    % accelerometer
    z_apri_acc_c = qRot(qRot(g,q_s_c), q_apri_c )/norm(qRot(qRot(g,q_s_c), q_apri_c ));
    % gyro
    z_apri_gyro_c = w_apri_c + gyro_bias_c;
    % Magnetometer
    z_apri_mag_c = qRot(qRot(mag_t,q_s_c), q_apri_c )/norm(qRot(qRot(mag_t,q_s_c), q_apri_c ));
    % Reaction Wheels
    z_apri_rw = rw_omega_apri;
    
    % Returning the apriori measurement vector
    z_apri = [z_apri_acc_c;
              z_apri_gyro_c;
              z_apri_mag_c;
              z_apri_rw];
    
    % Returning the apriori state vector
    x_apri = [q_apri_c;         % Includes q4
              w_apri_c;
              gyro_bias_c;
              rw_omega_apri];
    
    % Linearising and returning the system model and sensor model in x_apri
    [phi, H] = ekf_jacobian(x_apri, z_apri, Ji, J, q_s_c, q_s_tetra, dt);

    % propagating and returning the aposteriori error covariance matrix for
    % a apriori state error covariance matrix
    P_apri = phi * P_apo * phi' + Q;
end


% ======================================================================
%> @brief Fourth order Runge-Kutta
%>
%> @param J: Inertia matrix of the system
%> @param Ji: Inverse inertial matrix
%> @param w_apo_c:
%> @param act_torque: Actuator Torque
%> @param q_apo_c:
%> @param rw_omega_apo:
%> @param u: Input to actuator
%> @param dt:
%> @param rate:
%>
%> @retval x_out: Resulting state vector
% ======================================================================
function [x_out] = rk4(J, Ji, w_apo_c, act_torque, q_apo_c, rw_omega_apo, u, dt, rate)
    x = [q_apo_c;
         w_apo_c;
         rw_omega_apo];

    step_size = dt/rate;
    y = x;

    for i = 1:rate-1
        k1 = f(y, J, Ji, act_torque, u);
        k2 = f((y+0.5*step_size*k1), J, Ji, act_torque, u);
        k3 = f((y+0.5*step_size*k2), J, Ji, act_torque, u);
        k4 = f((y+step_size*k3), J, Ji, act_torque, u);

        y = y + (step_size*(k1 + 2*k2 + 2*k3 + k4))/6;
        y(1:4) = qunit(y(1:4));
    end
    x_out = y;
end


% ======================================================================
%> @brief System Model
%>
%> @param x: state vector
%> @param J: Inertia matrix of the system
%> @param Ji: Inverse inertial matrix
%> @param act_torque: Actuator Torque
%> @param u: Input to actuator
%>
%> @retval x_out: Resulting state vector
% ======================================================================
function [x_out] = f(x, J, Ji, act_torque, u)
    
    q_apo_c = x(1:4);
    w_apo_c = x(5:7);
    rw_omega_apo = x(8:11);
    
    % motor constants
    Kt = 1.81e-3;
    Ke = 1.81e-3;
    Rm = 4.44;
    bm = 6.3895e-08;
    Jrw = eye(4)*0.3811e-6;

    Lambda_m = eye(4,4) * ( (Kt * Ke) / Rm + bm );
    Gamma_m = eye(4,4)*(Kt/Rm);

    % system dynamics
    w_dot_c = Ji*(-skew_matrix(w_apo_c)*J*w_apo_c + act_torque);

    % system kinematics
    omega_w_apo = [skew_matrix(w_apo_c) w_apo_c;
                   -w_apo_c'             0];

    q_dot_c = 0.5*omega_w_apo*q_apo_c;

    % Motor dynamics
    rw_omega_dot = pinv(Jrw)*(-Lambda_m*rw_omega_apo + Gamma_m*u);

    x_out = [q_dot_c;
             w_dot_c;
             rw_omega_dot];
end


% ======================================================================
%> @brief Running the update step of the Multiplicative Extended Kalman Filter
%>
%> @param x_apri: Apriori estimated state vector
%> @param z_meas: Measurement vector
%> @param z_apri: Apriori estimated measurement vector
%> @param phi: Linearised system model, linear map
%> @param H: Linearised measurement linear map
%> @param P_apo_old: Previous aposteriori prediction error covariance matrix
%> @param Q: Proces noise covariance matrix
%> @param R: Measurement noise covariance matrix
%>
%> @retval x_apo: The apriori state estimate
%> @retval P_apo: Aporsteriori prediction error covariance matrix
% ======================================================================
function [x_apo, P_apo] = ekf_update(x_apri, z_meas, z_apri, H, P_apri, R, q_s_c)
    % Scale factor matrices for the sensor models
    RPM_RANGE = 16000;
    RPM2RAD_S = 2*pi/60;
    UNSIGNED_INT = 2^16;
    
    K_rpm = RPM_RANGE*RPM2RAD_S/UNSIGNED_INT; % + RPM uint16 resolution

    q_apri_full = x_apri(1:4);

    % Calculating the Kalman gain matrix
    S = pinv(H*P_apri*H' + R);
    K = P_apri*H'*S;

    % roatating z_meas into control frame
    z_meas_c = [qRot(z_meas(1:3), q_s_c); % accelerometer
                qRot(z_meas(4:6), q_s_c); % gyro
                qRot(z_meas(7:9), q_s_c); % magnetometer
                z_meas(10:13)*K_rpm];     % Reaction Wheels

    delta_x_apo = K*(z_meas_c - z_apri);
    delta_q_bar_apo = delta_x_apo(1:3);
    
    % [Humphreys,2002]
    tmp_dot = delta_q_bar_apo'*delta_q_bar_apo;
    if tmp_dot < 1
        delta_q_bar_apo_real = sqrt(1 - tmp_dot);
        q_update = [delta_q_bar_apo; delta_q_bar_apo_real];
    else
        delta_q_bar_apo_real = 1/sqrt(1 + tmp_dot);
        q_update = [delta_q_bar_apo*delta_q_bar_apo_real; delta_q_bar_apo_real];
        disp('tmp_dot way to large')
        disp(tmp_dot)
    end

    q_apo_t = qmult(q_apri_full,q_update);
    q_apo = qunit(q_apo_t);

    x_apo_tmp = x_apri(5:14) + delta_x_apo(4:13);
    
    x_apo = [q_apo;
             x_apo_tmp];

    % returning the aposteriori state vector error covariance matrix
    P_apo = ( eye(size(K*H, 1)) - K*H ) * P_apri;
end


% ======================================================================
%> @brief Performing a linearisation of the system model 
%>        and measurement model
%>
%> phi = d/dt * f(x,u) |(x=x_apri)
%> H = d/dt * h(x) |(x=x_apri)
%>
%> @param x_apri:
%> @param z_apri:
%> @param Ji:
%> @param J:
%> @param q_s_c:
%> @param q_s_tetra:
%> @param dt:
%>
%> @retval phi: Systen natrix linearised in x_apri
%> @retval H: Output matrix linearised in x_apri
% ======================================================================
function [phi,H] = ekf_jacobian(x_apri, z_apri, Ji, J, q_s_c, q_s_tetra, dt)
    q_apri_c = x_apri(1:4);
    w_apri_c = x_apri(5:7);
    gyro_bias_apri = x_apri(8:10);
    rw_omega = x_apri(11:14);

    S_w = skew_matrix(w_apri_c);
    S_Jw = skew_matrix(J*w_apri_c);

    % motor constants
    Kt = 1.81e-3;
    Ke = 1.81e-3;
    Rm = 4.44;
    bm = 6.3895e-08;
    Jrw = eye(4)*0.3811e-6;

    Lambda_m = eye(4,4) * ( ((Kt * Ke) / Rm) + bm );
    % Angles in the tetrahedron configuration
    alpha = deg2rad(60);
    beta = deg2rad(19.47);

    % tetrahedron projection matrix
    P_w_th = [cos(beta)    -cos(beta)*cos(alpha)   -cos(beta)*cos(alpha)    0;
              0             cos(beta)*cos(alpha/2) -cos(beta)*cos(alpha/2)  0;
             -sin(beta)    -sin(beta)              -sin(beta)               1];
         

    q_th_c = qmult(qinv(q_s_tetra),q_s_c);
    R_th_c = quat2rotm([q_th_c(4) q_th_c(1) q_th_c(2) q_th_c(3)]);
    
    rw_momentum = R_th_c*P_w_th*Jrw*rw_omega;
    skew_rw_momentum = skew_matrix(rw_momentum);
    
    f_w_q_partial = [-S_w        0.5*eye(3);
                     zeros(3)    Ji*((S_Jw - S_w*J) + skew_rw_momentum)];


    f_rw_partial = [                            zeros(3,4);
                    -Ji*(S_w*R_th_c*P_w_th*Jrw - R_th_c*P_w_th*Lambda_m)];
    
    rw_rw_partial = -pinv(Jrw)*Lambda_m;

    % The jacobian system matrix
    phi_t = [f_w_q_partial   zeros(6,3)   f_rw_partial;
             zeros(3,13);
             zeros(4,6)      zeros(4,3)   rw_rw_partial;];

    % ZOH discretisation
    phi = eye(size(phi_t,1)) + phi_t*dt;

    % the jacobian of the measurement matrix
    H = [2*skew_matrix(z_apri(1:3))  zeros(3)    zeros(3)  zeros(3,4);
         zeros(3)                    eye(3)      eye(3)    zeros(3,4);
         2*skew_matrix(z_apri(7:9))  zeros(3)    zeros(3)  zeros(3,4);
         zeros(4,3)                  zeros(4,3)  zeros(4,3)  eye(4)];
end


% ======================================================================
%> @brief Calculate the reaction wheel torque
%>
%> @param q_s_c:
%> @param q_s_tetra:
%> @param w_s_apo:
%> @param rw_omega_apo:
%> @param u:
%>
%> @retval tau_c_rw: three-dimensional torque vector in control reference frame
% ======================================================================
function tau_c_rw = reaction_wheels(q_s_c, q_s_tetra, w_s_apo, rw_omega_apo, u)
    % motor constants
    Kt = 1.81e-3;
    Ke = 1.81e-3;
    Rm = 4.44;
    bm = 6.3895e-08;
    Jrw = diag([0.4513 0.4513 0.3811]*1e-6);

    % Angles in the tetrahedron configuration
    alpha = deg2rad(60);
    beta = deg2rad(19.47);

    A1_q = [0.5773 -0.5773 0.4083 0.4083]';
    A2_q = [0.7887 0.2113 -0.1494 0.5577]';
    A3_q = [0.2113 0.7887 -0.5577 0.1494]';
    A4_q = [0 0 0.7071 0.7071]';

    omega_i_tetra = qRot(w_s_apo,q_s_tetra);
    q_th_c = qmult(qinv(q_s_tetra),q_s_c);

    omega_w1 = qRot(omega_i_tetra,qinv(A1_q)) + [0 0 rw_omega_apo(1)]';
    omega_w2 = qRot(omega_i_tetra,qinv(A2_q)) + [0 0 rw_omega_apo(2)]';
    omega_w3 = qRot(omega_i_tetra,qinv(A3_q)) + [0 0 rw_omega_apo(3)]';
    omega_w4 = qRot(omega_i_tetra,qinv(A4_q)) + [0 0 rw_omega_apo(4)]';

    tau_th_coi = (qRot(-cross(omega_w1,Jrw*omega_w1),A1_q)  +...
            qRot(-cross(omega_w2,Jrw*omega_w2),A2_q)) +...
            qRot(-cross(omega_w3,Jrw*omega_w3),A3_q)  +...
            qRot(-cross(omega_w4,Jrw*omega_w4),A4_q);

    % tetrahedron projection matrix
    P_w_th = [cos(beta)    -cos(beta)*cos(alpha)    -cos(beta)*cos(alpha)    0;
              0             cos(beta)*cos(alpha/2)  -cos(beta)*cos(alpha/2)  0;
             -sin(beta)    -sin(beta)               -sin(beta)               1];

    Lambda_m = eye(4)*(Kt*Ke/Rm + bm);

    Gamma_m = eye(4)*(Kt/Rm);

    tau_th_ctl = P_w_th*(-Lambda_m*rw_omega_apo + Gamma_m*u);

    tau_c_rw = -qRot((tau_th_coi + tau_th_ctl), q_th_c);
end
