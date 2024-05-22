function [res] = sun(t)
% SunV2 (SCT) wrapper: Calculate Sun position and distance
%
% Input:
%   t			Julian Date
%
% Output:
%   res(1,3)		ECI-based vector to Sun [m]



% Call function
[u,d]=SunV2(t);

% Place results in structure
res(1)=u(1)*d*1000;
res(2)=u(2)*d*1000;
res(3)=u(3)*d*1000;
