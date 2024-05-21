%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file contains the constants and other
% values used across the different designs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% operating point
h     = [2.23 2.23 2.23]*1e-3;       % momentum wheel bias
Sh    = [0 -h(3) h(2);h(3) 0 -h(1);-h(2) h(1) 0];
I     = diag([42.3 42.3 28.4]*1e-3); % inertia matrix for deployed situation