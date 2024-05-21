% Calculate angle between two grid index pairs.
% $Id: gridangle.m,v 1.1 2006/03/17 12:13:56 mnkr02 Exp $

function rho = gridangle(i1,j1,i2,j2,sy,sx);

[theta1,phi1] = idx2rad(i1,j1,sy,sx);
[theta2,phi2] = idx2rad(i2,j2,sy,sx);

rho = acos(sin(phi1)*sin(phi2)*cos(theta1-theta2)+cos(phi1)*cos(phi2));