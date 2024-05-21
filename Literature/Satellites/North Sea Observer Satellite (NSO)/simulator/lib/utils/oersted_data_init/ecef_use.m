function [sys,x0,str,ts] = ecef_use(t,x,u,flag)

% rot = ecef_use(theta, phi)
% Calculates rotation matrix rot(:,3,3) which transforms
% from an ECEF cartesian system to an USE (B_r, B_theta, B_phi)
% system. B-vector USE is rotated to B-vector ECEF.
% Input is B-vector, co-latitude theta(:) and longitude phi(:) in degrees

% Constants

switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
     [sys,x0,str,ts]=mdlInitializeSizes;
     
  %%%%%%%%%%%   
  % Outputs %
  %%%%%%%%%%%
  case 3,
    sys=mdlOutputs(t,x,u);
     
  %%%%%%%%%%%%%%%%%%%%%
  % Other valid flags %
  %%%%%%%%%%%%%%%%%%%%%
  case {1,2,4,9} 
    sys = []; 

  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%
  otherwise
    error(['Unhandled flag = ',num2str(flag)]);
end


%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%============================================================================
function [sys,x0,str,ts]=mdlInitializeSizes

sizes = simsizes;

sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0; 
sizes.NumOutputs     = 3; % ECEF B-vector
sizes.NumInputs      = 5; % USE B-vector, co-latitude and longitude
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1; % at least one sample time is needed

sys = simsizes(sizes); 


% initialize the initial conditions
x0  = [];

% str is always an empty matrix
str = [];

% initialize the array of sample time
ts  = [-1 0]; 

%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
function sys=mdlOutputs(t,x,u)

B(1) = u(1);
B(2) = u(2);
B(3) = u(3);
theta =  u(4);
phi =  u(5);

rad=pi/180;
ct = cos(theta*rad);
st = sin(theta*rad);
cp = cos(phi*rad);
sp = sin(phi*rad);

rot(1,1) = st.*cp;
rot(1,2) = +st.*sp;
rot(1,3) = +ct;
rot(2,1) = ct.*cp;
rot(2,2) = +ct.*sp;
rot(2,3) = -st;
rot(3,1) = -sp;
rot(3,2) = cp;
rot(3,3) = 0;

B(4) = rot(:,1)'*B(1:3)';
B(5) = rot(:,2)'*B(1:3)';
B(6) = rot(:,3)'*B(1:3)';

sys = [B(4:6)];
