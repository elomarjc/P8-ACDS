function [sys,x0,str,ts] = qestimater(t,x,u,flag,parm);
% Integrates a starting quaternion influenced by an angular velocity over
% time


  switch flag,
    %%%%%%%%%%%%%%%%
    case 0,
        [sys,x0,str,ts]=mdlInitializeSizes;
    %%%%%%%%%%%%%    
    %Derivatives%    
    %%%%%%%%%%%%%    
    case 1,
	 sys=mdlDerivatives(t,x,u);
    
     %%%%%%%%%    
    case 3,
        sys=mdlOutputs(t,x,u);
    %%%%%%%%%%%%%%%    
    %Ubrugte flags%  
    %%%%%%%%%%%%%%%    
    case {2,4,9},
    
    %%%%%%%%%%%%%%%%%%%
    %unexspected flags%    
    %%%%%%%%%%%%%%%%%%%    
    otherwise
        error(['Unhandled flag = ',num2str(flag)]);
end

function [sys,x0,str,ts]=mdlInitializeSizes

global firsttime;
firsttime = 0;

  sizes=simsizes;
  
  sizes.NumContStates = 4; % Quaternion state
  sizes.NumDiscStates = 0; 
  sizes.NumOutputs = 4; % Quaternion output
  sizes.NumInputs = 7; % w input and q initil condition
  sizes.DirFeedthrough = 7; 
  sizes.NumSampleTimes = 1; 

  sys = simsizes(sizes);
  

  x0 = [0 0 0 0];
  


  str = [];
  

  ts = [0 0];
  
% End of mdlInitialize



function sys=mdlDerivatives(t,x,u)

sys=0.5.*[x(4) -x(3) x(2) x(1);
         x(3) x(4) -x(1) x(2);
         -x(2) x(1) x(4) x(3);
         -x(1) -x(2) -x(3) x(4)]*[u(1); u(2); u(3); 0];

function sys=mdlOutputs(t,x,u);
global firsttime;
if(firsttime == 0)
  x(1) = u(4); 
  x(2) = u(5); 
  x(3) = u(6); 
  x(4) = u(7);
  firsttime = 1;
end

sys=x;