% PLOT_ALB_SPHERE Plot albedo data on a sphere.
%
% plot_alb_sphere(albedo)
%
% albedo is a 180x288 matrix.
%
% $Id: plot_alb_sphere.m,v 1.1 2006/03/17 12:13:57 mnkr02 Exp $

function h = plot_alb_sphere(albedo)

d2r = pi/180;

[sy sx] = size(albedo);
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
hp = surf(X,Y,Z,[albedo';albedo(:,1)']);
hold off;
shading('interp');
title('Albedo Data');
%colormap([[0:0.01:0.95]' [0:0.01:0.95]' [0:0.01:0.95]']);
view(127.5,16);

if nargout > 0
    % Return handle if requested
    h = hp;
end
