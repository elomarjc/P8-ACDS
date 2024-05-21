% CELLAREA_VERIFY Cell area verification. Calculates the Earth surface area
% by summing all grids and comparing to result from sphere surface from
% radius.
%
% cellarea_verify
%
% $Id: cellarea_verify.m,v 1.1 2006/03/17 12:13:56 mnkr02 Exp $

function ret = cellarea_verify;

global CONST;

area = zeros(180,288);

for i=1:180
	for j=1:288
		area(i,j) = cellarea(i,j,180,288);
	end
end

tarea = 4*pi*CONST.EMR^2;
aarea = sum(sum(area));
res = (tarea-aarea)/tarea;

fprintf('4*pi*r^2 = %.4g.\n',tarea);
fprintf('Accumulated area = %.4g\n',aarea);
fprintf('Accumulated area differs %2.2f%%.\n',res*100);

if nargout > 0
	ret = area;
end
