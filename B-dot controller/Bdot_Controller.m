%% setup
%magnetorquer_calculator;
%Stabilizing_Torque_calculator;
clc
%MagnetorquerPlant = Transferfunction;
dyear = decyear('20-may-2024','dd-mmm-yyyy');

BeField = load('BeField.mat');

%% orbital information
clc
n = 16470;
t = 0:1/n+1E4:1;
phi = 98*pi/180; %inclination
psi = sin(t); %longitude
theta = acos(sqrt(1- sin(psi).^2  *sin(phi)^2));  %lattitude
DataMagFieldBe = [psi;theta;zeros(size(t))];
Data = [psi;theta];
plot(t,Data)

%% world magnetic simulink block
clc
Data = out.Magfield;
%% Bdot data


DataBdotSatRot = out.DataBDotSatRot;
DataBdotMoment = out.DataBdotMoment;
minmaxBdotSatRot = out.MinMaxBdotSatRot;
minmaxBdotMoment = out.MinMaxBdotMoment;

%Rotation min/max lines
DataMaxRot = zeros(2.647E4,1);
DataMinRot = DataMaxRot;
DataMaxRot = DataMaxRot + 0.3*pi/180;
DataMinRot = DataMinRot - 0.3*pi/180;

%Moment min/max lines
DataMaxMom = zeros(2.647E4,1);
DataMinMom = DataMaxMom;
DataMaxMom = DataMaxMom + 105.9E-3;
DataMinMom = DataMinMom - 105.9E-3;


%% printing
figure(1);
plot(DataBdotSatRot), hold on, grid on
xline(n,'--'), 
plot(DataMaxRot), hold on
plot(DataMinRot),title("Satellite rotations"), xlabel("time [s]"), ylabel("angular velocity [rad/s]"), legend("X","Y","Z","req:T1.b","Max","Min", 'fontsize', 15), hold off

figure(2);
plot(DataBdotMoment), hold on, grid on
plot(DataMaxMom), hold on
plot(DataMinMom), title("Magnetic moment"), xlabel("time [s]"), ylabel("mag. moment [mA m^2]"), legend("X","Y","Z","Max","Min", 'fontsize', 15), hold off
%% test
