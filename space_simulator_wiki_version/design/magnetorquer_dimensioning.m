%**************************************************************************
% This file is used to calculate some parameters in the dimensioning of
% the magnetorquer, which are used in the magnetorquer modelling.
%
% Author: Group 06gr1032, modified by 09gr935
%**************************************************************************
clear all
clc

%%%Constants%%%
minB = 18e-6;                     %Tesla over equator 
                                  %http://www.magnetometer.org/mag-magnetic-field.php
Nreq=200e-9;                      %required torque from one coil

Duty_cycle2=0.88;                 %duty cycle on actuation. Coil must be off during measurements.

required_mm = ((Nreq/Duty_cycle2)/(minB)); %required magnetic moment.
voltage = 1.25;                   %maximum nominal votage of ADCS PWM
alpha_0 = 3.9e-3;                 %copper resistivity temp coef.
sigma_0 = 1.7e-8;                 %copper resistivity at 293K
ro = 8.92e3;                      %copper material density
mu_0 = 4*pi*10^-7;                %vacuum permeability

%%%Coil dimensions%%%
length = 75e-3;
width = 75e-3;

%%%Environment%%% -> values taken from ørsted data (solarpanel min/max temp)
temp_wc_i = 273-25;               %worst case for current at -25 deg C, as least resistance
temp = 273+85;                    %worst case for magnetic moment at 85 deg C, as most resistance;

%%%Variable parameters used to dimension the coil%%%
windings = 250;                   %coil windings
wire_dia = 0.13e-3;               %wire diameter
max_current = 0.040/2;            %maxcurrent divided by 2 for space qual. ELFA catalog 2006

%Resistivity
sigma = sigma_0*(1+alpha_0*(temp-293));     
sigma_wc_i = sigma_0*(1+alpha_0*(temp_wc_i-293));

%%%Coils%%%
circum = 2*length +2*width;                              %Coil circumference
R = (windings*circum*sigma)/((pi/4)*wire_dia^2);         %Resistance at 85 deg
I = voltage/R;                                           %Current at 85 deg
R_wc = (windings*circum*sigma_wc_i)/((pi/4)*wire_dia^2); %Resistance at -25 deg
I_wc = voltage/R_wc;                                     %Current at -25 deg
P = I_wc^2*R_wc*Duty_cycle2;                             %Worst case power consumption 
M = windings*circum*((pi/4)*wire_dia^2)*ro;              %Mass
A = length*width;                                        %Coil area
mm = windings*I*A*Duty_cycle2;                           %Magnetic moment
l = sqrt(A);                                             %Quadratic coil length eqv. 
L = (2*sqrt(2)*mu_0*windings^2*A)/(pi*l)                %Coil self inductance

%Printout%
disp(strcat('R=',num2str(R),' | I=',num2str(I),' | R_wc=',num2str(R_wc),' | I_wc=',num2str(I_wc)))
disp(strcat('P=',num2str(P),' | M=',num2str(M),' | mm=',num2str(mm),' | mm_req=',num2str(required_mm)))
mm_ok = ' PASSED';
if(mm < required_mm)                                     %Min. magnetic moment check
    mm_ok = strcat(' FAILED - Need aditional: ',num2str(required_mm-mm),' [Am^2]');
end
disp(strcat('Magnetic moment check = ',mm_ok))
i_ok = ' PASSED';
if(I_wc > max_current)                                   %Max current check
    v_max = R_wc*max_current;                            %Max allowed voltage
    R_add = voltage/max_current - R_wc;                  %Needed resistence
    i_ok = strcat(' FAILED - R_wc=',num2str(R_wc),...
        ' - Voltage limit=',num2str(v_max),' - Resistor needed=',num2str(R_add));
end
disp(strcat('Current check = ',i_ok))
Total_power=(6*0.005+6*P)*1000;
disp(strcat('Total worst case power consumption including drivers [mW] = ',num2str(Total_power))) %(all coils are saturated)
Total_mass=M*6*1000;
disp(strcat('Total mass of all coils [g]= ',num2str(Total_mass)))
disp(strcat('Dimensions [mm^2]= ',num2str(length*1000),'x ',num2str(width*1000)))
disp(strcat('Windings per coil= ',num2str(windings)))
disp(strcat('Coil diameter [mm]= ',num2str(wire_dia*1000)))

fill_factor=1.5; %the insulation takes up some space aswell (number taken from Lars)
Coil_cross_area=2*windings*(pi/4)*wire_dia^2*fill_factor;
disp(strcat('Coil cross sectional area [mm^2]= ',num2str(Coil_cross_area*1000*1000)))

disp(strcat('Max voltage pr. coil (controlled by PWM duty) [V]= ',num2str(voltage)))
disp(strcat('Duty cycle (controlled by enable pin) [%]= ',num2str(100*Duty_cycle2)))
distance=0.03; %e.g. distance from magnetorquer to magnetometer (Try not to exceed 100000nT)
Far_field=2*windings*mu_0*I_wc*length^2/(2*pi*(distance^2+length^2)/4*sqrt(distance^2+length^2/2)); % [Serway,p.965]
disp(strcat('Far field [nT]= ',num2str(10^9*Far_field),', at distance [mm]= ',num2str(1000*distance)))
tau_discharge=5*L/R_wc; %5 times time constant approx 99 discharged%
disp(strcat('Coil time constant (99% discharge time) [ms]= ',num2str(1000*tau_discharge)))
meas_time=(1-Duty_cycle2)*1/10-tau_discharge; %5 times time constant approx 99 discharged%
disp(strcat('Amount of time left for measurements [ms]= ',num2str(1000*meas_time)))

%% Coil with iron core %%

Nreq_ir=400e-9;         %required torque from one coil
Duty_cycle2_ir=0.27;    %duty cycle on actuation. Coil must be off during measurements.

required_mm_ir = ((Nreq_ir/Duty_cycle2_ir)/(minB)); %required magnetic moment.
voltage_ir = 1.25;      % maximum nominal votage of ADCS PWM
mu_r_ir = 1000;         % Relative permability is typically between 600-7600 and very temp dependent
%we should just find a material that has 1000 in permability minimum over
%whole temperature scale.
ro_ir = 7.874e3;        % Density of iron at 20 degress [kg/m^3]

%%%Iron core dimensions%%%
length_ir = 10e-3;
diameter_ir = 10e-3;

%%%Variable parameters used to dimension the coil%%%
windings_ir = 200;         %coil windings (it was tested that this is possible to have on 10x10mm core)
wire_dia_ir = 0.13e-3;     %wire diameter (with insulation 0.15e-3)
max_current_ir = 40e-3/2;  %maxcurrent divided by 2 for space qual. ELFA catalog 2006

%Resistivity
sigma_ir = sigma_0*(1+alpha_0*(temp-293));     
sigma_wc_i_ir = sigma_0*(1+alpha_0*(temp_wc_i-293));

%%%Coils%%%
circum_ir = (diameter_ir/2)^2*pi;                                   %Coil circumference
R_ir = (windings_ir*circum_ir*sigma)/((pi/4)*wire_dia_ir^2);        %Resistance at 85 deg
I_ir = voltage_ir/R_ir;                                             %Current at 85 deg
R_wc_ir = (windings_ir*circum_ir*sigma_wc_i)/((pi/4)*wire_dia_ir^2);%Resistance at -25 deg
I_wc_ir = voltage_ir/R_wc_ir;                                       %Current at -25 deg
P_ir = I_wc_ir^2*R_wc_ir*Duty_cycle2_ir;                            %Worst case power consumption
A_ir=(pi/4)*diameter_ir^2;
M_ir = windings_ir*circum_ir*((pi/4)*wire_dia_ir^2)*ro+((pi/4)*diameter_ir^2)*length_ir*ro_ir; %Mass
mm_ir = (mu_r_ir-1)*windings_ir*I_ir*((pi/4)*diameter_ir^2)*Duty_cycle2_ir; %Magnetic moment
L_ir =0;                                                            %Coil self inductance

%I_nes_ir=((Nreq_ir/Duty_cycle2_ir)/(minB));
I_nes_ir=(Nreq_ir/((mu_r_ir)*windings_ir*A_ir*minB));


%Printout%
disp('------------------------------------------------------------------------------------');
disp(strcat('R=',num2str(R_ir),' | I=',num2str(I_ir),' | R_wc=',num2str(R_wc_ir),' | I_wc=',num2str(I_wc_ir)))
disp(strcat('P=',num2str(P_ir),' | M=',num2str(M_ir),' | mm=',num2str(mm_ir),' | mm_req=',num2str(required_mm_ir)))
mm_ok = ' PASSED';
if(mm_ir < required_mm_ir)                                    %Min. magnetic moment check
    mm_ok = strcat(' FAILED - Need aditional: ',num2str(required_mm_ir-mm_ir),' [Am^2]');
end
disp(strcat('Magnetic moment check = ',mm_ok))
i_ok = ' PASSED';
if(I_wc_ir > max_current_ir)                                  %Max current check
    v_max = R_wc_ir*max_current_ir;                           %Max allowed voltage
    R_add = voltage_ir/max_current_ir - R_wc_ir;                 %Needed resistence
    I_mm = voltage_ir/(R_add+R_ir);
    i_ok = strcat(' FAILED - R_wc=',num2str(R_wc_ir),...
        ' - Voltage limit=',num2str(v_max),' - Resistor needed=',num2str(R_add));
end
disp(strcat('Current check = ',i_ok))
P_ir = max_current_ir^2*(R_ir+R_add)*Duty_cycle2_ir;
disp(strcat('Power consumption with R_add [mW] = ',num2str(P_ir*1e3))) %One coil
Total_power=(3*0.005+3*P_ir)*1000;
disp(strcat('Total worst case power consumption including drivers [mW] = ',num2str(Total_power))) %(all coils are saturated)
Total_mass=M_ir*3*1000+3*pi*length_ir*(diameter_ir/2)^2*ro_ir*1000;
disp(strcat('Total mass of all coils [g]= ',num2str(Total_mass)))
%disp(strcat('Dimensions [mm^2]= ',num2str(length_ir*1000),'x ',num2str(width*1000)))
disp(strcat('Windings per coil= ',num2str(windings_ir)))
%disp(strcat('Coil diameter [mm]= ',num2str(wire_dia*1000)))

%fill_factor=1.5; %the insulation takes up some space aswell (number taken from Lars)
%Coil_cross_area=2*windings*(pi/4)*wire_dia^2*fill_factor;
%disp(strcat('Coil cross sectional area [mm^2]= ',num2str(Coil_cross_area*1000*1000)))

disp(strcat('Max voltage pr. coil (controlled by PWM duty) [V]= ',num2str(voltage_ir)))
disp(strcat('Duty cycle (controlled by enable pin) [%]= ',num2str(100*Duty_cycle2_ir)))
distance=0.05; %e.g. distance from magnetorquer to magnetometer (Try not to exceed 100000nT)
%Far_field=2*windings*mu_0*I_wc*length^2/(2*pi*(distance^2+length^2)/4*sqrt(distance^2+length^2/2));
%disp(strcat('Far field [nT]= ',num2str(10^9*Far_field),', at distance [mm]= ',num2str(1000*distance)))
mm_ir2 = (mu_r_ir-1)*windings_ir*I_mm*((pi/4)*diameter_ir^2)*Duty_cycle2_ir; 
disp(strcat('mm_ir [Am^2]= ',num2str(mm_ir2)));
%disp(strcat('Try this max current [A]= ',num2str(I_nes_ir)));




