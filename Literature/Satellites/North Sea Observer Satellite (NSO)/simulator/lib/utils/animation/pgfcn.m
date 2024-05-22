function [sys,x0,str,ts] = TLEDecoder(t,x,u,flag);
  switch flag,
    % Initialize
    case 0,
        [sys,x0,str,ts]=mdlInitializeSizes;
    
    % Compute
    case 3,
        sys=mdlOutputs(t,x,u);

    % Unused
    case {1,2,4},
    
    % Simulation end
    case 9,
	 sys=mdlTerminate(t,x,u);

    % Faults
    otherwise
        error(['Unhandled flag = ',num2str(flag)]);
end


function [sys,x0,str,ts]=mdlInitializeSizes
  % Setup s-fcn size structure
  sizes=simsizes;
  sizes.NumContStates = 0;
  sizes.NumDiscStates = 3;
  sizes.NumOutputs = 1;
  sizes.NumInputs = 1;
  sizes.DirFeedthrough = 1;
  sizes.NumSampleTimes = 1;
  sys = simsizes(sizes);

  % Calculate total range
  startTime = str2num(get_param(bdroot,'Start time'));
  stopTime = str2num(get_param(bdroot,'Stop time'));
  pgRange = stopTime - startTime;

  % Prepare some nice text in the waitbar
  pgTxt = strcat(bdroot, ' running... Have patience...');
  pgTxt = strrep(pgTxt, '_', '-');

  % Initialize bar
  pgWaitbar = waitbar(0,pgTxt);

  % Store initialization results
  str = [];
  x0 = [pgWaitbar pgRange startTime];
  ts = [0 0];

function sys=mdlOutputs(t,x,u)
   % How far are we? u(1)=simtime, x(2)=range, x(3)=startOffset
   i = (u(1)-x(3)) / x(2);

   % Update waitbar with x(1)=handle
   waitbar(i, x(1));

   sys=i;

function sys=mdlTerminate(t,x,u)
   % Simply close our waitbar
   waitbar(1, x(1), 'Simulation finished');
   sys=[];

