%**************************************************************************
% This file plots the projected vectors in 3D, furthermore is also plots
% the vectors needed to calculate the Geometric Center of the area inclosed
% by the projected vectors.
% 
% The following .m files are needed for this file to run: arrow3.m, 
% nice3d.m and texts.m
% These files are included in the path: lib\utils\math_utils
%
% Author: Group 06gr1032
%**************************************************************************
format long;

nul = [0 0 0];
size = 12;      % Font size for the texts function.

% Projected vectors generated in the GCtest.mdl Simulink file.
p_x = px.signals.values(1,:);
p_y = py.signals.values(1,:);
p_z = pz.signals.values(1,:);
p_GC = pGC.signals.values(1,:);

% Calculating vectors used for off setting the projected vectors.
p_zx = p_x + p_z;
p_yx = p_y + p_x;
p_zxy = p_zx + p_y;

% Calculating the vector from [0 0 0] to the Geomatric Center.
vec_GC = (p_y + p_z + p_zx + p_yx + p_zxy)/6;

figure(1)
clf
% Plotting dimensions of the satellite also used as the reference frame.
arrow3(nul, [0.1 0 0], 'k'); nice3d % x-axis in meters.
hold on
texts([0.11 0 0], 'x',size)
arrow3(nul, [0 0.1 0], 'k'); nice3d % y-axis in meters.
texts([0 .11 0], 'y',size)
arrow3(nul, [0 0 .3], 'k'); nice3d  % z-axis in meters.
texts([0 0 .31], 'z',size)

% Plotting the projected vectors.
arrow3(nul, p_x, 'r'); nice3d
texts(p_x*1.2, 'px',size)
arrow3(nul, p_y, 'r'); nice3d
texts(p_y*1.2, 'py',size)
arrow3(nul, p_z, 'r'); nice3d
texts(p_z*1.2, 'pz',size)

% Plotting the offset vectors
arrow3(p_x, p_z, 'b'); nice3d
arrow3(p_z, p_x, 'b'); nice3d
arrow3(p_y, p_x, 'b'); nice3d
arrow3(p_x, p_y, 'b'); nice3d
arrow3(p_zx, p_y, 'b'); nice3d
arrow3(p_yx, p_z, 'b'); nice3d

% Plotting the GC projected from the reference frame.
arrow3(nul, p_GC, 'g'); nice3d

% Plotting the last three vectors used for calculating the Geometric Center
% of the area inclosed by the projected vectors. Besides these three the
% vectors p_y and p_z are used in the calculation.
arrow3(nul, p_yx, 'y'); nice3d
arrow3(nul, p_zx, 'y'); nice3d
arrow3(nul, p_zxy, 'y'); nice3d

% Plotting the calculated GC.
arrow3(nul, vec_GC, 'm'); nice3d

axis off

% It is noticed that the calculated and projected GC coinside which means
% that the projected GC is also the center of pressure of the exposed area.