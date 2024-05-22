% DOUBLE_PLOT_EPR Plot epr data on subplot with dual view.
%
% double_plot_epr(epr)
%
% epr is either a epr struct or epr matrix.
%
% $Id: double_plot_epr.m,v 1.1 2006/03/17 12:13:50 mnkr02 Exp $

function double_plot_epr(epr);

h = figure;
subplot(1,2,1);
plot_epr(epr);
subplot(1,2,2);
plot_epr(epr);
view(90,0);

return