
if exist('RunFile','var')==0
    clc
    close all
    clear
    WCDTorque
end

%% Init
if exist('RunFile','var')==1

%http://www.powermagnetshop.de/pd-1250611506.htm?categoryId=3
Bi=1.21; % Intrinsic induction or remanence [T]
dpmag=2e-3; % Diameter of magnet [m]
hpmag=1e-3; % Heigth of magnet [m]
Vpmag=pi*(dpmag/2)^2*hpmag; % Volume of magnet [m^3]

%mmag=0.005; % [A*m^2]
%nmag=mmag*60000e-9
%fsat=(1/(2*pi))*sqrt(nmag/Isat2(3,3))
%Periode=1/fsat
%Vpm=((mmag*mu0)/Bi)*1e9

mpmag=(Bi*Vpmag)/mu0 % Magnetic moment [Am^2]

Bffpmag=((mu0*mpmag)/(2*(0.05^2+dpmag^2)^(3/2)))*1e9 % [nT]

npmag_min=norm(mpmag)*norm(minB) % Torque [Nm]
npmag_max=norm(mpmag)*norm(maxB) % Torque [Nm]

fsat_min=(1/(2*pi))*sqrt(npmag_min./Isat) % Oscillation frequency [Hz]
fsat_max=(1/(2*pi))*sqrt(npmag_max./Isat) % Oscillation frequency [Hz]

Periode=1./fsat_min % Oscillation periode [s]
Periode=1./fsat_max % Oscillation periode [s]

%orbits=(Isat2(3,3)*((InitDetumble)/nmag))/Torbit

end

