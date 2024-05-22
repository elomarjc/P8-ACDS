% ALBEDO_ALTITUDE Calculate subsolar albedo between two altitudes (above
% earth) at specified position.
%
% result = albedo_altitude(az,pa,alt1,alt2,epr,n,type);
%
% az and pa is the azimuth and polar angle of the position in ECEF frame.
% Albedo is calculated from alt1 to alt2 with n intermediate calculations
% (defaults to 10). Note that altitudes must be given in meters. Type
% specifies the result type. '%' for percent values otherswise result is
% W/m^2. Use 'p' for plot and '%p' for percentage plots. The first 5
% parameters are required.
%
% $Id: albedo_altitude.m,v 1.1 2006/03/17 12:13:49 mnkr02 Exp $

function result = albedo_altitude(az,pa,alt1,alt2,epr,n,type);

CONST.EMR = 6371.01e3;
CONST.AU = 149598000000;
CONST.AM0 = 1366.9;

if nargin < 5
  error('Not enough input parameters');
end

if nargin < 6
  type = '';
	n = 10;
end

if nargin < 7
  if ischar(n)
    type = n;
    n = 10;
  else
    type = '';
  end
end

maxalt = max(alt1,alt2);
minalt = min(alt1,alt2);
step = floor((maxalt-minalt)/n);
result = zeros(n+1,1);
i = 1;
h = waitbar(0,'Calculating Albedos...');
for alt = minalt:step:maxalt
	a = albedo([az;pa;CONST.EMR+alt],[az;pa;CONST.AU],epr);
	result(i) = sum(sum(a));
	waitbar(i/(n+1),h);
	i = i + 1;
	drawnow;
end

if ~isempty(findstr('%',type))
	result = result./CONST.AM0;
	yname = 'Albedo [E.U.]';
else
	yname = 'Albedo [W/m^2]';
end

close(h);

if ~isempty(findstr('p',type))
  plot([minalt:step:maxalt]./1000,result);
  title('Subsolar Albedo');
  xlabel('Altitude [km]');
  ylabel(yname);
end
