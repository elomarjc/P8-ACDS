% PLOT_EPR_SPHERE Plot epr data on a sphere.
%
% plot_epr_sphere(epr_parm)
%
% epr_parm is a struct containing a 180x288 matrix in a
% field named data or a 180x288 matrix.
%
% $Id: plot_epr_sphere.m,v 1.1 2006/03/17 12:13:57 mnkr02 Exp $

function h = plot_epr_sphere(epr_parm)

if ~isfield(epr_parm,'data')
  epr.data = epr_parm;
else
  epr = epr_parm;
end

d2r = pi/180;

[sy sx] = size(epr.data);
dy = 180/sy;
dx = 360/sx;

phi = [-90+dy/2:dy:90-dy/2]'.*d2r;
theta = [-180+dx/2:dx:180-dx/2]'.*d2r;

X = cos(theta)*cos(phi)';
Y = sin(theta)*cos(phi)';
Z = ones(sx+1,1)*sin(phi)';

X = [X;X(1,:)];
Y = [Y;Y(1,:)];

hp = surf(X,Y,Z,zeros(sx+1,sy));
hold on;
hp = surf(X,Y,Z,[epr.data';epr.data(:,1)']);
hold off;
shading('interp');
title('TOMS Reflectivity Data');
%colormap([[0:0.01:0.95]' [0:0.01:0.95]' [0:0.01:0.95]']);
view(127.5,16);

if nargout > 0
    % Return handle if requested
    h = hp;
end
