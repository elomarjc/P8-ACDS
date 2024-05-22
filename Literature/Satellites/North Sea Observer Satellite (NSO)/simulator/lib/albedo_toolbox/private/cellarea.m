% CELLAREA Calculate area of TOMS cell. Input is epr matrix
% indices.
%
% A = cellarea(i,j,sy,sx)
%
% i,j are the indices of the syxsx epr data matrix.
%
% $Id: cellarea.m,v 1.1 2006/03/17 12:13:56 mnkr02 Exp $

function A = cellarea(i,j,sy,sx);

CONST.EMR = 6371.01e3;
CONST.d2r = pi/180;

[theta,phi] = idx2rad(i,j,sy,sx);

% Grid size
dphi = (180/sy)*CONST.d2r;
dtheta = (360/sx)*CONST.d2r;

% Diagonal points
phimax = phi + dphi/2;
phimin = phi - dphi/2;

A = CONST.EMR^2*dtheta*(cos(phimin)-cos(phimax));
