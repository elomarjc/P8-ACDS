% Transform TOMS EPR matrix indices to radians
% $Id: idx2rad.m,v 1.1 2006/03/17 12:13:56 mnkr02 Exp $

function [theta,phi] = idx2rad(i,j,sy,sx);

CONST.d2r = pi/180;

dx = 2*pi/sx;
dy = pi/sy;

phi = pi-dy/2-(i-1)*dy;
theta = (j-1)*dx-pi+dx/2;
