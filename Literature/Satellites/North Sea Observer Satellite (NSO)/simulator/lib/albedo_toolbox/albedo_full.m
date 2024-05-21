% ALBEDO_FULL Calculate albedo array for 100% sunlight (zenith) at all
% satellite positions at a given altitude. The resulting albedo is the
% maximum albedo for all satellite positions.
%
% a = albedo_full(altitude,epr);
%
% altitude is the altitude of the satellite in meters. epr is the
% reflectivity data.
%
% $Id: albedo_full.m,v 1.1 2006/03/17 12:13:50 mnkr02 Exp $

function a = albedo_full(altitude,epr);

CONST.EMR = 6371.01e3;
CONST.AU = 149598000000;

h = waitbar(0,'Calculating full albedo matrix...');
for i = 1:180
	for j = 1:288
		[sat_theta sat_phi] = idx2rad(i,j);
		satsph = [sat_theta,sat_phi,CONST.EMR+altitude];
		sunsph = [sat_theta,sat_phi,CONST.AU];
		a(i,j) = sum(sum(albedo(satsph,sunsph,epr)));
		waitbar(i*j/(180*288),h);
	end
end

close(h);
