% EARTHFOV Field of view on earth by spherical coordinates.
%
% result = earthfov(satsph, epr, type)
%
% satsph is the vector to the satellite in ECEF frame and spherical
% coordinates. epr is the reflectivity data used for plotting the fov.
% type is 'p' for 3D plot and 's' for spherical. If epr or type are
% unspecified no plot is generated. epr is not specified, epr must be a
% vector of the latitudal (epr(1)) and longitudal (epr(2)) resolution of
% the output data.
%
% $Id: earthfov.m,v 1.1 2006/03/17 12:13:50 mnkr02 Exp $

function result = earthfov(satsph,epr,type);

CONST.EMR = 6371.01e3;
CONST.d2r = pi/180;

% LEO shortcut
if satsph(3) < CONST.EMR
	satsph(3) = satsph(3) + CONST.EMR;
end

% Circle value
OUTVALUE = 1;

% Discretization parameters
if isfield(epr,'data')
  [sy sx] = size(epr.data);
else 
  sy = epr(1);
  sx = epr(2);
end
dx = 2*pi/sx;
dy = pi/sy;
result = zeros(sy,sx);

% Small Circle Center
theta0 = satsph(1);
phi0 = satsph(2);

% FOV on earth
rho = acos(CONST.EMR/satsph(3));

for i = 1:sy
	for j = 1:sx        
		[theta, phi] = idx2rad(i,j,sy,sx);
		% Radial distance
		rd = acos(sin(phi0)*sin(phi)*cos(theta0-theta)+cos(phi0)*cos(phi));
		if rd <= rho
			result(i,j) = OUTVALUE;
		end
	end
end

if nargin > 1 && isfield(epr,'data')
  if nargin > 2
    plot_epr(mask(epr.data,result),type);
  else
    plot_epr(mask(epr.data,result));
  end
end
