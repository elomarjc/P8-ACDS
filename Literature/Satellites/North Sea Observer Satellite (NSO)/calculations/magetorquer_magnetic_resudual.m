%**************************************************************************
% This file is used to calculate some parameters in the dimensioning of
% the magnetorquer, which are used in the magnetorquer modelling.
%
% Author: Group 06gr1032
%**************************************************************************
clear all
clc

%%%Constants%%%
required_mm = 130.94e-3; %required magnetic moment.
voltage = 5;             %maximum nominal votage of ADCS PWM
alpha_0 = 3.9e-3;        %copper resistivity temp coef.
sigma_0 = 1.7e-8;        %copper resistivity at 293K
ro = 8.92e3;             %copper material density
mu_0 = 4*pi*10^-7;       %vacuum permeability

%%%Coil dimensions%%%
xy_length = 210e-3;
xy_width = 70e-3;
z_length = 80e-3;
z_width = 80e-3;

%%%Environment%%%
temp_wc_i = 273-25;     %worst case for current at -25 deg C, as least resistance
temp = 273+70;          %worst case for magnetic moment at 70 deg C, as most resistance;

%%%Variable parameters used to dimension the coil%%%
windings_xy = 435;      %coil windings
windings_z = 1;       %coil windings
wire_dia = 0.1e-3;      %wire diameter
max_current = 93e-3/2;  %maxcurrent divided by 2 for space qual. ELFA catalog 2006

%Resistivity
sigma = sigma_0*(1+alpha_0*(temp-293));     
sigma_wc_i = sigma_0*(1+alpha_0*(temp_wc_i-293));

%%%Coils on x and y panels%%%
circum_xy = 2*xy_length +2*xy_width;                              %Coil circumference
R_xy = (windings_xy*circum_xy*sigma)/((pi/4)*wire_dia^2);         %Resistance at 70 deg
I_xy = voltage/R_xy;                                              %Current at 70 deg
R_wc_xy = (windings_xy*circum_xy*sigma_wc_i)/((pi/4)*wire_dia^2); %Resistance at -25 deg
I_wc_xy = voltage/R_wc_xy;                                        %Current at -25 deg
P_xy = I_wc_xy^2*R_xy;                                            %Worst case power consumption 
M_xy = windings_xy*circum_xy*((pi/4)*wire_dia^2)*ro;              %Mass
A_xy = xy_length*xy_width;                                        %Coil area
mm_xy = windings_xy*I_xy*A_xy;                                    %Magnetic moment
l_xy = sqrt(A_xy);                                                %Quadratic coil length eqv. 
L_xy = (2*sqrt(2)*mu_0*windings_xy^2*A_xy)/(pi*l_xy);             %Coil self inductance

%Printout%
disp(strcat('R_xy=',num2str(R_xy),' | I_xy=',num2str(I_xy),' | I_wc_xy=',num2str(I_wc_xy)))
disp(strcat('P_xy=',num2str(P_xy),' | M_xy=',num2str(M_xy),' | mm_xy=',num2str(mm_xy)))
disp(strcat('L_xy=',num2str(L_xy)))
mm_ok = ' PASSED';
if(mm_xy < required_mm)                                    %Min. magnetic moment check
    mm_ok = ' FAILED';
end
disp(strcat('Magnetic moment check = ',mm_ok))
i_ok = ' PASSED';
if(I_wc_xy > max_current)                                  %Max current check
    v_max_xy = R_wc_xy*max_current;                        %Max allowed voltage
    R_add_xy = voltage/max_current - R_wc_xy;              %Needed resistence
    i_ok = strcat(' FAILED - R_wc_xy=',num2str(R_wc_xy),...
        ' - Voltage limit=',num2str(v_max_xy),' - Resistor needed=',num2str(R_add_xy));
end
disp(strcat('Current check = ',i_ok))


disp('....................................................') %Following is eqv to above

%%%z coil%%%
circum_z = 2*z_length + 2*z_width;;
R_z = (windings_z*circum_z*sigma)/((pi/4)*wire_dia^2);
I_z = voltage/R_z;
R_wc_z = (windings_z*circum_z*sigma_wc_i)/((pi/4)*wire_dia^2);
I_wc_z = voltage/R_wc_z;
P_z = I_wc_z^2*R_z;
M_z = windings_z*circum_z*((pi/4)*wire_dia^2)*ro;
A_z = z_length*z_width;
mm_z = windings_z*I_z*A_z;
l_z = sqrt(A_z); 
L_z = (2*sqrt(2)*mu_0*windings_z^2*A_z)/(pi*l_z);

%Printout%
disp(strcat('R_z=',num2str(R_z),' | I_z=',num2str(I_z),' | I_wc_z=',num2str(I_wc_z)))
disp(strcat('P_z=',num2str(P_z),' | M_z=',num2str(M_z),' | mm_z=',num2str(mm_z)))
disp(strcat('L_z=',num2str(L_z)))
mm_ok = ' PASSED';
if(mm_z < required_mm)
    mm_ok = ' FAILED';
end
disp(strcat('Magnetic moment check = ',mm_ok))
i_ok = ' PASSED';
if(I_wc_z > max_current)
    v_max_z = R_wc_z*max_current;
    R_add_z = voltage/max_current - R_wc_z;
    i_ok = strcat(' FAILED - R_wc_z=',num2str(R_wc_z),...
        ' - Voltage limit=',num2str(v_max_z),' - Resistor needed=',num2str(R_add_z));
end
disp(strcat('Current check = ',i_ok))