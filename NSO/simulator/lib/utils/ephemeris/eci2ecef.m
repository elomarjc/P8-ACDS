function Qeci2ecef = eci2ecef(t)
% Calculate ECI to ECEF rotation quartenion by Julian Date.
%
% q = eci2ecef(t)

% Calculates the ECI->ECEF (z-axis) quartenion given a Julian date.
% Output is within 0 - 2pi boundary.
% Originally based on eci2ecef.m by raf@control.auc.dk.

% Rotation of the Earth with respect to the mean equinox is referred to as
% mean sidereal time (ST). The Earth rotates 360 deg in one sidereal day
% which is 23h 56m 04.09053s [almanac, b6] in mean solar time. 
% The resulting Earth precession rate is therefore 2*pi/23h 56m 04.09053s
% Apart from inherent motion of the equinox, due to precession and nutation,
% ST is a direct measure of the diurnal rotation of the Earth. As a fixpoint
% for the Earth rotation the Greenwich median transit of the equinox at
% 0 Jan 1997 is used (JD2450448.5). According to [almanac, b15] the
% transit takes place at 17h 18m 21.8256s (JD0.7211).
% The Earth Precession is denoted "omega"
% The sidereal transit is denoted "s"

% Seconds on a day
daysec = 86400;

% Earth precession [rad/sec]
omega = 7.292115855377075e-05;

% Sidereal epoch (raf's)
s = 6.23018356e+04/daysec + 2450448.5;

% Sidereal epoch (jkpe01's)
%s = 2453096.22108;

% Time for a revolution
T_r = 2*pi/omega;

% Number of revolutions since epoch
N_r = round( ((t-s)*daysec) / T_r);

% Time into a revolution
T_n = (t-s)*daysec - N_r*T_r;

% Calculate angle with 2pi-boundary
psi = mod( omega*T_n, 2*pi );

% Represent result as a quartenion
Qeci2ecef = [ 0;
	      0;
	      sin(psi/2);
	      cos(psi/2) ];
