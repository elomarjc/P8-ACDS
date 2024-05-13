% MEAN_REFLECTIVITY Calculate mean Earth reflectivity of epr data.
% Each data point is weighed with respect to the cell area of the
% associated cell.
%
% total_refl = mean_reflectivity(epr);
%
% epr is an epr struct.
%
% $Id: mean_reflectivity.m,v 1.1 2006/03/17 12:13:51 mnkr02 Exp $

function total_refl = mean_reflectivity(epr);

CONST.EMR = 6371.01e3;

total_area = 4*pi*CONST.EMR^2;
total_refl = 0;
for i=1:180
	for j=1:288
		weight = cellarea(i,j)/total_area;
		total_refl = total_refl + weight*epr.data(i,j);
	end
end
