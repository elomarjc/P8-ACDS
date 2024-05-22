% Used to define the liniar state space system that is used in the
% verification of the liniear model.

%Constants
I=diag([42.3 42.3 28.4]/1e3); % inertia matrix for deployed situation
w=zeros(3);
Sw=[0 -w(3) w(2);w(3) 0 -w(1);-w(2) w(1) 0];
Iw=I*w;
SIw=[0 -Iw(3) Iw(2);Iw(3) 0 -Iw(1);-Iw(2) Iw(1) 0];
h=zeros(1,3);
Sh=[0 -h(3) h(2);h(3) 0 -h(1);-h(2) h(1) 0];

%Defining the linear model
A_model_ver=[-Sw 0.5*eye(3) zeros(3)
              zeros(3) inv(I)*(SIw-Sw*I+Sh) -inv(I)*Sw 
              zeros(3) zeros(3) zeros(3)];
B_model_ver=[zeros(3) zeros(3)
              inv(I) -inv(I)
              zeros(3) eye(3)];
C_model_ver=eye(9);
D_model_ver=zeros(9,6);


%Max torque on one axis
maxDist = 1599.1e-9; %Nm
maxMagMoment = 202.99e-3; %Am^2
maxBField = 70e-6; %T - http://www.magnetometer.org/mag-magnetic-field.php
maxMW = 0.817e-3; %Nm from datasheet. Friction neglected.
%------
maxTorque = maxDist+maxMW+maxMagMoment*maxBField;