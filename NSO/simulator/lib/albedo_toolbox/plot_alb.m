% PLOT_ALB Plot albedo data.
%
% plot_alp(albedo,type)
%
% albedo is a 180x288 matrix. Type is 'p' for 3D plot or 's' for spherical. 
%
% $Id: plot_alb.m,v 1.1 2006/03/17 12:13:51 mnkr02 Exp $

function h = plot_alb(albedo,type);

if nargin > 1 && strcmp(type,'s')
	plot_alb_sphere(albedo);
	return;
end

[sy sx] = size(albedo);
dy = 180/sy;
dx = 360/sx;

lat = [-90+dy/2:dy:90-dy/2]';
lon = [-180+dx/2:dx:180-dx/2]';

hp = surf(lon,lat,albedo);
view(0,90);
axis([-180+dx 180-dx -90+dy 90-dy]);
shading('interp');
title('Albedo Data');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
zlabel('Albedo [W/m^2]');
%colormap([[0:0.01:0.95]' [0:0.01:0.95]' [0:0.01:0.95]']);

% return handles, if requested
if nargout > 0
    h = hp;
end
