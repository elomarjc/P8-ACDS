% ALBEDO_CSUN Calculate albedo array for constant sunlight (instantanious
% albedo) for all satellite positions at given altitude.
%
% a = albedo_csun(altitude,sunsph,epr,type);
%
% altitude is the altitude of the albedo calculations in meters. sunsph is the
% position of the Sun in ECEF spherical coordinates. epr is the
% reflectivity data. Type is 'p' for 3D plot and 's' for spherical. If
% unspecified no plot is generated.
%
% $Id: albedo_csun.m,v 1.1 2006/03/17 12:13:50 mnkr02 Exp $

function a = albedo_csun(altitude,sunsph,epr,type);

CONST.EMR = 6371.01e3;

if nargin > 3
  % Sunlit elements
  vis = earthfov(sunsph);
  plot_epr(mask(epr.data,vis),type);
end

h = waitbar(0,'Calculating instantanious constant sun albedo matrix...');
for i = 1:180
	for j = 1:288
		[sat_theta sat_phi] = idx2rad(i,j);
		satsph = [sat_theta,sat_phi,CONST.EMR+altitude];
		a(i,j) = sum(sum(albedo(satsph,sunsph,epr)));
		waitbar(i*j/(180*288),h);
  end
end

close(h);
