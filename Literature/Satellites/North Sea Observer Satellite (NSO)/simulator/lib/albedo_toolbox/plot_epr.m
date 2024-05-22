% PLOT_EPR Plot epr data.
%
% plot_epr(epr_parm,type)
%
% epr_parm is a struct containing a 180x288 matrix in a
% field named data or a 180x288 matrix. Type is 'p' for 3D plot or 's' for
% spherical.
%
% $Id: plot_epr.m,v 1.1 2006/03/17 12:13:52 mnkr02 Exp $

function h = plot_epr(epr_parm,type);

if ~isfield(epr_parm,'data')
  epr.data = epr_parm;
else
  epr = epr_parm;
end
  
if nargin > 1 && strcmp(type,'s')
	plot_epr_sphere(epr);
	return;
end

[sy sx] = size(epr.data);
dy = 180/sy;
dx = 360/sx;

lat = [-90+dy/2:dy:90-dy/2]';
lon = [-180+dx/2:dx:180-dx/2]';

hp = surf(lon,lat,epr.data);
view(0,90);
axis([-180+dx 180-dx -90+dy 90-dy]);
shading('interp');
title('TOMS Reflectivity Data');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
zlabel('Reflectivity [E.U.]');
%colormap([[0:0.01:0.95]' [0:0.01:0.95]' [0:0.01:0.95]']);

% return handles, if requested
if nargout > 0
    h = hp;
end
