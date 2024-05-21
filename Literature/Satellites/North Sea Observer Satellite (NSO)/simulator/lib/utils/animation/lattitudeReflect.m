function [sys,x0,str,ts] = lattitudeReflect(t,x,u,flag);
  switch flag,
    % Initialize
    case 0,
        [sys,x0,str,ts]=mdlInitializeSizes;
    
    % Compute
    case 3,
        sys=mdlOutputs(t,x,u);

    % Unused
    case {1,2,4, 9},
    
    % Faults
    otherwise
        error(['Unhandled flag = ',num2str(flag)]);
end


function [sys,x0,str,ts]=mdlInitializeSizes
  % Setup s-fcn size structure
  sizes=simsizes;
  sizes.NumContStates = 0;
  sizes.NumOutputs = 2;
  sizes.NumInputs = 2;
  sizes.DirFeedthrough = 1;
  sizes.NumSampleTimes = 1;

  % Ugly hack, but we actually load the everything into the discrete
  % state stuff
  load albedotab;
  len = length(ads_c_mean);

  % Insert length and store data in state
  x0(1) = len;
  for i=1:len
      x0(1+i) = ads_c_mean(i);
      x0(1+i+len) = ads_c_devi(i);
  end;

  % Store sizes struct
  sizes.NumDiscStates = 2*len + 1;
  sys = simsizes(sizes);

  % Store initialization results
  str = [];
  ts = [0 0];

function sys=mdlOutputs(t,x,u)
   % Inputs:
   % u(1) - lattitude of SC
   % u(2) - field of view

   % Lattitudes per index in data
   lpi = 180/x(1);

   % Lattitude and FOV in degrees
   la = 180/pi*u(1);
   Cfov = 180/pi*u(2);

   % Start with zero albedo and deviation
   Rref = 0;
   Rdev = 0;
   
   % Number of steps to take, and the starting fov lattitude
   steps = abs(2*floor(Cfov / lpi));
   fov = -1*Cfov;

   % Step from -fov to fov
   for n = 0:steps
   
      % Get albedo gain at current point. The function will
      % calculate the correct index including 0/180 wrap-hack
      [rt dt] = albedoL(la + fov, x);

      % FOV area compensation for both reflectivity and deviation
      compen = cos(fov/Cfov * pi/2);

      if compen > 1
	 error('compen > 1')
      end
      if rt > 1
	 error('rt > 1')
      end

      Rref = Rref + rt*compen/(steps+1);
      Rdev = Rdev + dt*compen/(steps+1);

      % Increment fov by the resolution
      fov = fov + lpi;
   end

   % Just store result, and we're done
   sys = [ Rref Rdev ];

function [ref,dev] = albedoL(la,x)
  
   % Lattitudes per index in data
   lpi = 180/x(1);

   % Because of the area compensation addition, we'll have to wrap
   % around both 0 and 180 degrees
   if la < 0
      la = abs(la);
   end
   if la > 180
      la = 180 - rem(la, 180);
   end

   % Calculate index of lattitude, x(1) is the no. of lattitudes stored
   idx = floor(abs(la/lpi)) + 2;

   % Return reflectivity and standard deviation
   ref = x(idx);
   dev = x(idx + x(1));
