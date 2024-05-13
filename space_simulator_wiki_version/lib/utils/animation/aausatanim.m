function [sys,x0,str,ts] = aausatanim(t,x,u,flag,Config)
% Based on
% SAEROANIM S-Function for displaying 6DoF trajectories
% Copyright 1990-2002 The MathWorks, Inc.
% J.Hodgson  
% $Revision: 1.2 $  $Date: 2006/05/29 16:22:02 $

switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
     [sys,x0,str,ts]=mdlInitializeSizes(Config);

  %%%%%%%%%%%%%%%
  % Derivatives %
  %%%%%%%%%%%%%%%
  case {1 , 3, 9},
     sys=[];
     
  %%%%%%%%%%
  % Update %
  %%%%%%%%%%
  case 2,
     sys = [];
     if Config.Animenable
        mdlUpdate(t,x,u,Config);
     end

  %%%%%%%%%%%%%%%%%%%%%%%
  % GetTimeOfNextVarHit %
  %%%%%%%%%%%%%%%%%%%%%%%
  case 4,
    sys=mdlGetTimeOfNextVarHit(t,x,u,Config);

 
otherwise
  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%

   error(['Unhandled flag = ',num2str(flag)]);

end

% end sanim
%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts]=mdlInitializeSizes(Config)

%
% Set Sizes Matrix
%
sizes = simsizes;

sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 0;
sizes.NumInputs      = 14;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

%
% initialise the initial conditions
%
x0  = [];

%
% str is always an empty matrix
%
str = [];

%
% initialise the array of sample times
%
ts  = [-2 0]; % variable sample time

if ~Config.Animenable
   return
end

%
% Initialise Figure Window
%

   h_f=findobj('type','figure','Tag','aausatiianim');
   
   if isempty(h_f)
     h_anim=figure;
   else
     h_anim=h_f;
   end

   set(h_anim,'name','AAUSAT-II Animation', ...
           'renderermode','auto','resize','on', ...
           'position',[150 150 700 700],'clipping','off', ...
           'Tag','aausatiianim','renderer','zbuffer');
   
   if ~isempty(h_anim)
      h_del = findobj(h_anim,'type','axes');
      delete(h_del);
      figure(h_anim);
   end
%
% Initialize Axes
%
   handle.axes(1)=axes;
   axis(Config.axes);
   set(handle.axes(1),'visible','on','xtick',[],'ytick',[],'ztick',[],'box','off', ...
           'dataaspectratio',[1 1 1], ...
           'projection','pers', ...
           'units','normal', ...
           'position',[0.1 0.1 0.75 0.75], ...
           'Color',[0 0 0], ...
           'drawmode','fast', ...
           'clipping','off');
   handle.light = light('Position',[0 2e14 0],'color',[1 1 0.5]);
%
% Initialize Trajectory 
%
    if(Config.Trajectory) 
        handle.line(1) = line(0,0,0);
        set(handle.line(1),'linestyle','-','color',[1 0.7 0.3],'erasemode','nor','userdata',0,'clipping','off');
    end
      
%
% Draw in Central body Position
% Karl's magic sphere generator
%   
    res=(2*pi)/Config.PlanetSlices;
    prd=Config.PlanetSlices/2+1;
    vrt=0;
    face=0;
    i=0;
    for long=0:res:2*pi-res
        i=i+1;
        for lat=pi/2:-res:-pi/2
            vrt=vrt+1;
            rxy=Config.radius*cos(lat);
            zz=Config.radius*sin(lat);
            xx=rxy*cos(long);
            yy=rxy*sin(long);
            verts(vrt,:)=[xx yy zz];
            if (mod(vrt,prd)~=0)
                face=face+1;
                if (i<Config.PlanetSlices)
                    fce(face,:)=[vrt vrt+1 vrt+1+prd vrt+prd];
                else
                    fce(face,:)=[vrt vrt+1 mod(vrt,prd)+1 mod(vrt,prd)];
                end
            end
        end
    end
	P.Vert = verts;
	P.faces = fce;

    % Now, map that Earth image onto the shpere, yeaahhhhh!!!! :-)
    if(Config.UseMap)
        im = imread(Config.PlanetMapFile);
        im = im2double(im);
        i=0;
        j=0;
        tcolors=ones(Config.PlanetSlices*(Config.PlanetSlices/2),3);
        for i=1:1:Config.PlanetSlices
            for j=1:1:Config.PlanetSlices/2
                tcolors(((i-1)*Config.PlanetSlices/2)+j,:) = im(j,i,:);
            end
        end
        handle.target = patch('vertices',P.Vert,'faces',P.faces,'FaceVertexCData',tcolors,'FaceColor','flat');
        set(handle.target,'erasemode','nor','edgecolor','none','clipping','off','userdata',[0]);
        rotate(handle.target,[0 0 1],Config.PreRotation);
    else
        handle.target = patch('vertices',P.Vert,'faces',P.faces);
        set(handle.target,'erasemode','nor','edgecolor',[0 0.3 1],'facecolor',[0 0.3 1],'clipping','off','userdata',[0]);
    end
    alpha(1)

    if(Config.ShowECEF)
        axislines = [0 1 1-0.1 1-0.1 1
                     0 0 0     0     0
                     0 0 0.05 -0.05  0];
        axislength = 1.3*Config.radius;
        
        ex = axislines*axislength;
        ey = [axislines(2,:);axislines(1,:);axislines(3,:)]*axislength;
        ez = [axislines(3,:);axislines(2,:);axislines(1,:)]*axislength;
            
        handle.efxaxis = line(ex(1,:), ex(2,:), ex(3,:));
        handle.efyaxis = line(ey(1,:), ey(2,:), ey(3,:));
        handle.efzaxis = line(ez(1,:), ez(2,:), ez(3,:));
	
        set(handle.efxaxis,'color',[0 1 1]);
        set(handle.efyaxis,'color',[1 0 1]);
        set(handle.efzaxis,'color',[1 1 0]);
    end

    if(Config.ShowInertial)
        axislines = [0 1 1-0.1 1-0.1 1
                     0 0 0     0     0
                     0 0 0.05 -0.05  0];
        axislength = 1.5*Config.radius;       
	
        ix = axislines*axislength;
        iy = [axislines(2,:);axislines(1,:);axislines(3,:)]*axislength;
        iz = [axislines(3,:);axislines(2,:);axislines(1,:)]*axislength;
            
        handle.xaxis = line(ix(1,:), ix(2,:), ix(3,:));
        handle.yaxis = line(iy(1,:), iy(2,:), iy(3,:));
        handle.zaxis = line(iz(1,:), iz(2,:), iz(3,:));
	
        set(handle.xaxis,'color',[0 1 1]);
        set(handle.yaxis,'color',[1 0 1]);
        set(handle.zaxis,'color',[1 1 0]);
        
        text(axislength,0.05*Config.radius,-0.05*Config.radius,'X','FontSize',12,'color',[1 1 1]);
        text(0.05*Config.radius,axislength,-0.05*Config.radius,'Y','FontSize',12,'color',[1 1 1]);
        text(0.05*Config.radius,-0.05*Config.radius,axislength,'Z','FontSize',12,'color',[1 1 1]);
    end
               
%
% Draw in Satellite
%
    if(Config.UseSCMap)
        cur = 1;
        fce = 1;
        face = ones((Config.CraftSlices+1)^2,2);
        one = ones((Config.CraftSlices+1)^2,1);
        squares = ones(Config.CraftSlices^2,4);
        P.tcolor = ones(Config.CraftSlices^2*6,3);
        for i=-Config.CraftSlices:2:Config.CraftSlices
            for j=-Config.CraftSlices:2:Config.CraftSlices
                face(cur,:)=[i j];
                if((j<Config.CraftSlices)&&(i<Config.CraftSlices))
                    squares(fce,:) = [cur cur+1 cur+Config.CraftSlices+2 cur+Config.CraftSlices+1];
                    fce=fce+1;
                end
                cur=cur+1;
            end
        end
        face=face/Config.CraftSlices;
        P.sides = [one       face(:,2) -3*face(:,1)
                 face(:,1) one       -3*face(:,2)
                 face(:,1) face(:,2) 3*one
                 -one      face(:,2) -3*face(:,1)
                 face(:,1) -one      -3*face(:,2)
                 face(:,1) face(:,2) -3*one      ];
        P.faces = [squares(:,:)
                 squares(:,:)+ones(Config.CraftSlices^2,4).*(Config.CraftSlices+1)^2
                 squares(:,:)+ones(Config.CraftSlices^2,4).*2*(Config.CraftSlices+1)^2
                 squares(:,:)+ones(Config.CraftSlices^2,4).*3*(Config.CraftSlices+1)^2
                 squares(:,:)+ones(Config.CraftSlices^2,4).*4*(Config.CraftSlices+1)^2
                 squares(:,:)+ones(Config.CraftSlices^2,4).*5*(Config.CraftSlices+1)^2];
        P.sides = P.sides.*Config.craft/2;
        P.tcolor = rand(Config.CraftSlices^2*6,3);
        im = imread('sat15.jpg');
        im = im2double(im);
        i=0;
        j=0;
        for i=1:1:Config.CraftSlices*6
            for j=1:1:Config.CraftSlices
                P.tcolor((i-1)*Config.CraftSlices+j,:) = im(i,j,:);
            end
        end
	
        handle.craft = patch('vertices',P.sides,'faces',P.faces,'FaceVertexCData',P.tcolor,'FaceColor','flat');
        set(handle.axes(1),'userdata',P.sides);
    else
        P.Vert = Config.craft/2*[-1 -1 -3;1 -1 -3;1 1 -3;-1 1 -3;-1 -1 3;1 -1 3;1 1 3;-1 1 3];
        P.faces = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];
        P.tcolor = [1 0 0;0 1 1;1 0 1;0 1 0;0 0 1;1 1 0];
        handle.craft = patch('vertices',P.Vert,'faces',P.faces,'FaceVertexCData',P.tcolor,'FaceColor','flat');
        set(handle.axes(1),'userdata',P.Vert);
    end
    set(handle.craft,'erasemode','nor','edgecolor','none','clipping','off');

%
% Space craft coordinate system   
%
    if(Config.ShowSCFrame)
        axislines = [0 1 1-0.1 1-0.1 1
                     0 0 0     0     0
                     0 0 0.05 -0.05  0];
        axislength = 1.5*Config.craft;       

        handle.cx = axislines*axislength;
        handle.cy = [axislines(2,:);axislines(1,:);axislines(3,:)]*axislength;
        handle.cz = [axislines(3,:);axislines(2,:);axislines(1,:)]*axislength*1.5;
        
        handle.cxaxis = line(handle.cx(1,:), handle.cx(2,:), handle.cx(3,:));
        handle.cyaxis = line(handle.cy(1,:), handle.cy(2,:), handle.cy(3,:));
        handle.czaxis = line(handle.cz(1,:), handle.cz(2,:), handle.cz(3,:));
	
        set(handle.cxaxis,'color',[0 1 1]);
        set(handle.cyaxis,'color',[1 0 1]);
        set(handle.czaxis,'color',[1 1 0]);
    end
%
% Orbit coordinate system
%

    if(Config.ShowOrbitFrame)
        handle.oxaxis = line([0 0],[0 0],[0 0]);
        handle.oyaxis = line([0 0],[0 0],[0 0]);
        handle.ozaxis = line([0 0],[0 0],[0 0]);
	
        set(handle.oxaxis,'color',[0 1 1]);
        set(handle.oyaxis,'color',[1 0 1]);
        set(handle.ozaxis,'color',[1 1 0]);
    end
    

%
% Set Handles of graphics in Figure UserData
%   
   set(h_anim,'userdata',handle);

%
%=============================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=============================================================================
%
function mdlUpdate(t,x,u,Config,count)

%
% Obtain Handles of Figure Objects
%
    handle = get(findobj('type','figure','Tag','aausatiianim'),'userdata');

    if isempty(findobj('type','figure','Tag','aausatiianim'))
     %figure has been manually closed
     return
    end

%
% Form Transformation Matrix

	q0 = u(4);
	q1 = u(5);
	q2 = u(6);
	q3 = u(7);
           
    attitude = [1-2*(q1^2+q2^2) 2*(q0*q1-q2*q3)  2*(q0*q2+q1*q3)
                2*(q0*q1+q2*q3) 1-2*(q0^2+q2^2)  2*(q1*q2-q0*q3)
                2*(q0*q2-q1*q3) 2*(q1*q2+q0*q3)  1-2*(q0^2+q1^2)];
%
% Update Craft Object 
%
	vert = get(handle.axes(1),'userdata');
	[a,b]=size(vert);
	dum =attitude*vert'+u(1:3,1)*ones(1,a);
	set(handle.craft,'vertices',dum');
	
    % Rotate SC coodinatesystem
	if(Config.ShowSCFrame)
		cx = attitude*handle.cx+u(1:3,1)*ones(1,5);
		cy = attitude*handle.cy+u(1:3,1)*ones(1,5);
		cz = attitude*handle.cz+u(1:3,1)*ones(1,5);
		set(handle.cxaxis,'XData',cx(1,:),'YData',cx(2,:),'ZData',cx(3,:));
		set(handle.cyaxis,'XData',cy(1,:),'YData',cy(2,:),'ZData',cy(3,:));
		set(handle.czaxis,'XData',cz(1,:),'YData',cz(2,:),'ZData',cz(3,:));
    end
	
    % Make orbit coordinate system
    if(Config.ShowOrbitFrame)
        norm = sqrt(u(1)^2+u(2)^2+u(3)^2);
		xx = u(1)/norm;
		xy = u(2)/norm;
		xz = u(3)/norm;
    	ox = -[0 xx
               0 xy
               0 xz];
		
		yx = u(2)*u(10)-u(3)*u(9);
		yy = u(3)*u(8)-u(1)*u(10);
		yz = u(1)*u(9)-u(2)*u(8);
		norm = sqrt(yx^2+yy^2+yz^2);
		yx = yx/norm;
		yy = yy/norm;
		yz = yz/norm;
    	oy = [0 yx
              0 yy
              0 yz];
		
		zx = u(2)*yz-u(3)*yy;
		zy = u(3)*yx-u(1)*yz;
		zz = u(1)*yy-u(2)*yx;
		norm = sqrt(zx^2+zy^2+zz^2);
		zx = zx/norm;
		zy = zy/norm;
		zz = zz/norm;
    	oz = -[0 zx
               0 zy
               0 zz];
           
        basis = [ox(:,2) oy(:,2) oz(:,2)];
        axislines = [0 1 1-0.1 1-0.1 1
                     0 0 0     0     0
                     0 0 0.05 -0.05  0];
        axislength = 1.5*Config.craft;
        rx = basis*axislines*axislength+u(1:3,1)*ones(1,5);
        ry = basis*[axislines(2,:);axislines(1,:);axislines(3,:)]*axislength+u(1:3,1)*ones(1,5);
        rz = basis*[axislines(3,:);axislines(2,:);axislines(1,:)]*axislength+u(1:3,1)*ones(1,5);
        
        set(handle.oxaxis,'XData',rx(1,:),'YData',rx(2,:),'ZData',rx(3,:));
		set(handle.oyaxis,'XData',ry(1,:),'YData',ry(2,:),'ZData',ry(3,:));
		set(handle.ozaxis,'XData',rz(1,:),'YData',rz(2,:),'ZData',rz(3,:));
    end
	
%
% Rotate the Earth
%
    angle = u(14)*180/pi;
    oldangle = get(handle.target,'userdata');
    rotate(handle.target,[0 0 1],angle-oldangle);
    set(handle.target,'userdata',angle);
    % Rotate ECEF coodinatesystem
	if(Config.ShowECEF)
        rotate(handle.efxaxis,[0 0 1],angle-oldangle);
        rotate(handle.efyaxis,[0 0 1],angle-oldangle);
        rotate(handle.efzaxis,[0 0 1],angle-oldangle);
    end

%
% Update Line Objects
%
    if(Config.Trajectory)
        x1 = get(handle.line(1),'Xdata');  
        x2 = get(handle.line(1),'Ydata');
        x3 = get(handle.line(1),'Zdata');
        init = get(handle.line(1),'userdata');
        if init
           x1 = [x1 u(1)];
           x2 = [x2 u(2)];
           x3 = [x3 u(3)];
        else
           x1 = [u(1)];
           x2 = [u(2)];
           x3 = [u(3)];
           set(handle.line(1),'userdata',1);
        end
        set(handle.line(1),'Xdata',x1,'Ydata',x2,'Zdata',x3);	
    end                   
%
% Set position of target view to Target
%
switch Config.camera_view 

case 1,			% Fixed Observer Position
   set(handle.axes(1),'cameraupvector',[0 0 1], ...
      'cameraposition',Config.camera_pos, ...
      'cameratarget',[0 0 0], ...
      'cameraviewangle',Config.view);
   
case 2,			% Relative Position View
   set(handle.axes(1),'cameraupvector',[0 0 -1], ...
           'cameraposition',u(1:3).*10+Config.camera_pos'.*0.2, ...
           'cameratarget',u(1:3), ...
           'cameraviewangle',Config.view);
end
set(handle.light, 'Position',[u(11) u(12) u(13)]);

%
% Force MATLAB to Update Drawing
%
   drawnow

   if(Config.MakeMovie)
     frame = getframe(gcf);
     P = frame2im(frame);
     directory = 'test/animation/movie/';
     number = 1000+num2str(ceil((t+1)/Config.update));
     extension = '.jpg';
     imwrite(P,strcat(directory,number,extension),'jpeg');
    end
   
   
% end mdlUpdate

%
%=============================================================================
% mdlGetTimeOfNextVarHit
% Return the time of the next hit for this block.  
%=============================================================================
%
function sys=mdlGetTimeOfNextVarHit(t,x,u,Config)
    
    sys = ceil(t/Config.update)*Config.update+Config.update;

% end mdlGetTimeOfNextVarHit

