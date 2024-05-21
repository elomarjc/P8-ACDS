% interface to msis86 routines
function [sys,x0,str,ts] = Qcontinuum(t,x,u,flag)

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
sizes.NumOutputs     = 1; % density
sizes.NumInputs      = 2; % 
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

time = u(1);
day = u(2);

DAY = day - 10; %oersted data start: _10._ february

seconds = DAY*86400 + time;

sys = seconds;


