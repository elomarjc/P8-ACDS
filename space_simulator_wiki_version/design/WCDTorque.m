
%% Check, Close and Clean


clear
clc
close all

%% Set falg
RunFile=1;

%% Parameters

G=6.672e-11;
mE=5.9742e24;
SatLength=100e-3;
SatHeigth=110e-3;
SatWidth=100e-3;
BoardLength=87e-3;
BoardWidth=87e-3;
ERadiusMean=6371e3;

CD=2.0;
Fsolar=1366;
c=2.99792458e8;
Crk=1.5;
nSolar=1;
nBoard=0; % Not valid using grund and Vcc plane
ISolarCell=487e-3*2;
SatStackPower=2.5;
SatStackVoltage=3.3;
muE=7.82579e15; % Earth's magnetic dipole streng calculated using [Wertz1]
mu0=(4*pi)*1e-7;
norbits=3;
InitDetumble=pi/18; % 10 [degree/s]
GoM=[0.05,0.05,0.055]';
mu=G*mE;
Prmf=Fsolar/c;

maxB=48000e-9; % 600 [km] orbit
minB=18000e-9; % 600 [km] orbit
AltitudeSat=600e3;
AeroAltitude=[100; 110; 120; 130; 140; 150; 160; 180; 200; 250; 300; 350; 
              400; 450; 500; 600; 700; 800; 900; 1000];
AeroRho=[5.297e-7; 9.661e-8; 2.438e-8; 8.484e-9; 3.845e-9; 2.070e-9; 1.244e-9; 
         5.464e-10; 2.789e-10; 7.248e-11; 2.418e-11; 9.158e-12; 3.725e-12; 
         1.585e-12; 6.967e-13; 1.454e-13; 3.614e-14; 1.170e-14; 5.245e-15; 3.019e-15];
rho=AeroRho(16); %1.454e-13


Torbit=2*pi*sqrt((ERadiusMean+AltitudeSat)^3/mu); % Orbit time [s]
V=(2*pi*(ERadiusMean+AltitudeSat))/Torbit; % Velocity in heigth of 700 [km]
EndDetumble=(4*pi)/Torbit;
omegaPoint=atan(V/AltitudeSat); 
%% Open files
SatelliteInertia

%plot altitude dependent torques
for k=15:20
    
    rho=AeroRho(k);
    AeroDragMax
    torque_drag(k-14)=m_torque;

    AltitudeSat=1000*AeroAltitude(k);
    GravityTorqueMax
    torque_gravity(k-14)=Ngg;  
    
    RadiationTorqueMax
    torque_radiation(k-14)=m_torque_r;
    
end
figure
hold on
hej=polyfit(AeroAltitude(15:20),torque_drag(:),5);
t=[500:1:900];
plot(t,(hej(1)*t.^5+hej(2)*t.^4+hej(3)*t.^3+hej(4)*t.^2+hej(5)*t+hej(6))*10^(9),'--k','LineWidth',2)
plot(AeroAltitude(15:19),torque_gravity(1:5)*10^(9),'-.k','LineWidth',2);
plot(AeroAltitude(15:19),torque_radiation(1:5)*10^(9),'-k','LineWidth',2);
fontsize=12;
xlabel('Altitude [km]','FontSize',fontsize);
ylabel('Torque [nNm]','FontSize',fontsize);
legend('Aerodynamic torque','Gravity gradient torque','Solar radiation torque');


%find worst case torques for 600 km height
rho=AeroRho(16);
AltitudeSat=600e3;
GravityTorqueMax
AeroDragMax
RadiationTorqueMax
MagneticResidualTorque
PermanentMagnetTorque
TumblePointTorque


disp('Total Disturbance Torque (High Field Stregnth)[Nm]')
nDistTotal=(Ngg+m_torque+max(m_torque_r)+Nmagsc_max)

disp('Total Torque With Pointing (High Field Strength) [Nm]')
nTotal=nDistTotal+max([Npoint,Ndet])+npmag_max

disp('Total Disturbance Torque (Low Field Stregnth)[Nm]')
nDistTotal=(Ngg+m_torque+max(m_torque_r)+Nmagsc_min)

disp('Total Torque With Pointing (Low Field Strength) [Nm]')
nTotal=nDistTotal+max([Npoint,Ndet])+npmag_min


%% Clear flag
clear RunFile