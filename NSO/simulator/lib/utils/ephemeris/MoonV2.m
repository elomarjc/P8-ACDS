function [u, r] = MoonV2( jD, rSc )

%-------------------------------------------------------------------------------
%   Generate the moon vector in the earth-centered inertial frame
%   or, if a spacecraft vector is input, in the spacecraft
%   centered frame. This is the moderate precision model. 
%-------------------------------------------------------------------------------
%   Form:
%   [u, r] = MoonV2( jD, rSc )
%-------------------------------------------------------------------------------
%
%   Inputs
%   ------
%   jd        (:)     Julian date
%   rsc       (:)     Spacecraft vector in the ECI frame (km)
%
%   -------
%   Outputs
%   -------
%   u         (3,:)   Unit moon vector
%   r         (:)     Distance from origin to moon (km) 
%
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
%	References: Montenbruck, O., T.Pfleger, Astronomy on the Personal
%                 Computer, Springer-Verlag, Berlin, 1991, pp. 103-111.
%-------------------------------------------------------------------------------
%	 Copyright 1993 Princeton Satellite Systems, Inc. All rights reserved.
%-------------------------------------------------------------------------------

if( nargin < 1 )
  jD = Date2JD;
end

% Days from J2000.0

T    = (jD - 2451545.0)/36525;
tpi  = 2*pi;

% The mean arguments neglecting long-period oscillations
% which have a maximum magnitude of 88.7"

L0   =     rem(0.60643382+(1336.85522467-3.130e-6*T).*T,1);
L    = tpi*rem(0.37489701+(1325.55240982+2.565e-5*T).*T,1);
LS   = tpi*rem(0.99312619+(  99.99735956-4.400e-7*T).*T,1);
D    = tpi*rem(0.82736186+(1236.85308708-3.970e-6*T).*T,1);
F    = tpi*rem(0.25909118+(1342.22782980-8.920e-6*T).*T,1);

TD   = 2*D;
TF   = 2*F;
TL   = 2*L;
LPLS = L+LS;

STF  = 412*sin(TF);
SLS  =     sin(LS);

% In arcseconds

DL   = 22640*sin(L)-4586*sin(L-TD)+2370*sin(TD)+769*sin(TL)-668*SLS-STF -212*sin(TL-TD)-206*sin(LPLS-TD)+192*sin(L+TD)-165*sin(LS-TD)-125*sin(D)-110*sin(LPLS)+148*sin(L-LS)-55*sin(TF-TD);

% Conversion from arcseconds to radians

k    = pi/(180*3600);

S    = F + (DL+STF+541*SLS)*k;
H    = F - TD;
N    = -526*sin(H) + 44*sin(L+H) - 31*sin(H-L) - 23*sin(LS+H) + 11*sin(H-LS) - 25*sin(F-TL) + 21*sin(F-L);

lam  = tpi*rem(L0+DL/1296e3,1);
beta = (18520*sin(S) + N)*k;

cb   = cos (beta);
sb   = sin (beta);
cl   = cos (lam);
sl   = sin (lam);

u = [cb.*cl;0.91748*cb.*sl-0.39778*sb;0.39778*cb.*sl+0.91748*sb];

% Account for parallax

if ( nargout == 2 | nargin == 2 ),
  dp = 186.5398*cos(L)+34.3117*cos(L-TD)+28.2333*cos(TD)+10.1657*cos(TD) -0.3997*cos(LS)-0.0124*cos(TF);
  sp = 0.999953253*(3422.7+dp)*k;
  r  = 6378.137./sp;
end


% Account for parallax
%---------------------
if ( nargin == 2 ),
  u = [r.*u(1,:) - rSc(1,:);...
       r.*u(2,:) - rSc(2,:);...
       r.*u(3,:) - rSc(3,:)];
  r = Mag(u);
  u = Unit(u);
end

  
