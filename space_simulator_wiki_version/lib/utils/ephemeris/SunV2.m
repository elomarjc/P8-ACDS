function [u, r] = SunV2( jD, rSc )

%-------------------------------------------------------------------------------
%   Generate the sun vector in the earth-centered inertial frame. Will
%   output the distance to the sun from the earth if two output
%   arguments are given.
%-------------------------------------------------------------------------------
%   Form:
%   [u, r] = SunV2( jD, rSc )
%-------------------------------------------------------------------------------
%
%   ------
%   Inputs
%   ------
%   jD        (:)     Julian date
%   rSc       (3,:)   Spacecraft vector in the ECI frame (km)
%
%   -------
%   Outputs
%   -------
%   u         (3,:)   Unit sun vector
%   r         (:)     Distance from origin to sun (km)
%
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
%   References: Montenbruck, O., T.Pfleger, Astronomy on the Personal
%               Computer, Springer-Verlag, Berlin, 1991, p. 36.
%-------------------------------------------------------------------------------
%   Copyright 1993 Princeton Satellite Systems, Inc. All rights reserved.
%-------------------------------------------------------------------------------

if( nargin < 1 )
  jD = [];
end

% Today
%------
if( isempty(jD) )
  jD = Date2JD;
end

twoPi = 2*pi;

% Days from J2000.0
%------------------
T     = (jD - 2451545.0)/36525;

% Mean anomaly
%-------------
g     = rem(0.993133 + 99.997361*T,1);
dLam  = 6893*sin( twoPi*g ) + 72*sin( twoPi*2*g );
lam   = twoPi*rem(0.7859453 + g + (6191.2*T + dLam)/1296e3,1);

% Obliquity of ecliptic
%----------------------
obOfE   = 0.408982 - 2.269e-04*T;

sLam  = sin(lam);

u     = [cos(lam); cos(obOfE).*sLam; sin(obOfE).*sLam];

if ( nargin == 2 | nargout == 2 ),
  r = (1.0014 - 0.01671*cos(g*pi/180) - 0.00014*cos(2*g*pi/180))*149600e3;
end

% Account for parallax
%---------------------
if ( nargin == 2 )
  u = [r.*u(1,:) - rSc(1,:);...
       r.*u(2,:) - rSc(2,:);...
       r.*u(3,:) - rSc(3,:)];
  r = Mag(u);
  u = Unit(u);
end
