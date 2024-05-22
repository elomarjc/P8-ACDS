%**************************************************************************
% This file is used to plot some illustrations used in the magnetorquer 
% modelling.
%
% The following .m files are needed for this file to run: arrow3.m, nice3d.m
% and texts.m
% These files are included in the path: ./
%
% Author: Group 06gr1032
%**************************************************************************
close all;
clear all;
clc;
format compact;

%%Prints 0 = off, 1 = on
printEPS = 1; %print to eps
forces = 1;
actuationForces = 0;
if(actuationForces)
    forces = 0;
end
outputDir = '';
outputName = 'magnetorquer_illu';
if(forces)
    outputName = 'magnetorquer_forces';
elseif(actuationForces)
    outputName = 'magnetorquer_actuation_forces';
end

I = 5 % 50 mA;
B = [.06 .06 .06]*0.7;

x = [0.02 0 0];
y = [0 0.02 0];
z = [0 0 0.02];

nul = [.05 -.05 0];
bnul = [0 0 0];

% Coil vectors
Sx1 = [-.10 0 0];
Sy1 = [0 .10 0];
Sx2 = [.10 0 0];
Sy2 = [0 -.10 0];

Sy1_start = nul;
Sx1_start = Sy1_start + Sy1;
Sy2_start = Sx1_start + Sx1;
Sx2_start = Sy2_start + Sy2;

% Force vectors
ISx1xB = cross(I*Sx1, B);
ISy1xB = cross(I*Sy1, B);
ISx2xB = cross(I*Sx2, B);
ISy2xB = cross(I*Sy2, B);

% Plot torquer and B-field vector
figure
hold on
arrow3(Sy1_start,Sy1)
if(forces==0 && actuationForces==0)
texts(Sy1_start + Sy1/4  + [0 -.008 -.005], 'S_{\alpha_1}', 12)
end
arrow3(Sx1_start,Sx1)
if(forces==0  && actuationForces==0)
texts(Sx1_start + Sx1/4 + [-.008  0 .005] , 'S_{\beta_1}', 12)
end
arrow3(Sy2_start,Sy2)
if(forces==0 && actuationForces==0)
texts(Sy2_start + Sy2/4 + [0 -.008 -.005] , 'S_{\alpha_2}', 12)
end
arrow3(Sx2_start, Sx2)
if(forces==0 && actuationForces==0)
texts(Sx2_start + Sx2/4 + [-.008 0 -.01], 'S_{\beta_2}', 12)
end
arrow3(bnul,B, 'r')
texts(B*1.01 + bnul, 'B', 12)

% Plot coordinate system vectors

% arrow3([0 0 0], x, 'k')
% texts([0 0 0] + x*1.3, 'X', 12)
% arrow3([0 0 0], y, 'k')
% texts([0 0 0] + y*1.3, 'Y', 12)
% arrow3([0 0 0], z, 'k')
% texts([0 0 0] + z*1.3, 'Z', 12)


% Plot force vectors
if(forces)
arrow3(Sx1_start + Sx1/2, ISx1xB, 'g')
texts(Sx1_start + 3*Sx1/4 + [-.000  0 .005] , 'S_{\beta_1}', 12)
texts(Sx1_start + Sx1/2 + ISx1xB*1.3 + [0.001 -0.005 0.009], 'F_{\beta_1}', 12)

arrow3(Sy1_start + Sy1/2,ISy1xB, 'g')
texts(Sy1_start + 3*Sy1/4  + [0 -.000 -.006], 'S_{\alpha_1}', 12)
texts(Sy1_start + Sy1/2 + ISy1xB*1.3 + [0.001 -0.015 0.015], 'F_{\alpha_1}', 12)

arrow3(Sx2_start + Sx2/2, ISx2xB, 'g')
texts(Sx2_start + 3*Sx2/4 + [-.008 0.006 -.01], 'S_{\beta_2}', 12)
texts(Sx2_start + Sx2/2 + ISx2xB*1.3 + [0.001 0.015 -0.01], 'F_{\beta_2}', 12)

arrow3(Sy2_start + Sy2/2,ISy2xB, 'g')
texts(Sy2_start + 3*Sy2/4 + [-.006 -.008 .008] , 'S_{\alpha_2}', 12)
texts(Sy2_start + Sy2/2 + ISy2xB*1.3 + [0.001 0.015 -0.01], 'F_{\alpha_2}', 12)
end

% % B projected to the yz plane
% By = B.*[0 1 1]
% 
% % B projected to the xz plane
% Bx = B.*[1 0 1]
% 
% theta_x = acos(dot(By,y)/(norm(By)*norm(y)))
% theta_y = acos(dot(Bx,x)/(norm(Bx)*norm(x)))
% deg_x = theta_x/(2*pi)*360
% deg_y = theta_y/(2*pi)*360
% 
% % Torque around X-axis
% F_x1 = norm(ISx1xB)*cos(pi+theta_x)*[0 0 1]
% F_y1 = norm(ISy1xB)*cos(pi+theta_y)*[0 0 1]
% F_x2 = norm(ISx2xB)*cos(theta_x)*[0 0 1]
% F_y2 = norm(ISy2xB)*cos(theta_y)*[0 0 1]
% 
% % Transformation matrix that maps to the plane spanned by x and z
% Axz=[x' z'];
% Qxz=Axz*(Axz'*Axz)^(-1)*Axz'
% 
% % Transformation matrix that maps to the plane spanned by y and z
% Ayz=[y' z'];
% Qyz=Ayz*(Ayz'*Ayz)^(-1)*Ayz'

% Fq_x1 = Qxz*ISx1xB'
% Fq_y1 = Qyz*ISy1xB'
% Fq_x2 = Qxz*ISx2xB'
% Fq_y2 = Qyz*ISy2xB'

% Reflection matrix that maps to the plane with normal vector x
X_norm = x'/norm(x);
Qr_yz = eye(3)-X_norm*X_norm'


% Reflection matrix that maps to the plane with normal vector y
Y_norm = y'/norm(y);
Qr_xz = eye(3)-Y_norm*Y_norm';

F_y1 = Qr_yz*ISy1xB';
F_x1 = Qr_xz*ISx1xB';
F_y2 = Qr_yz*ISy2xB';
F_x2 = Qr_xz*ISx2xB';

% Plot actuation force vectors
if(actuationForces)
arrow3(Sx1_start + Sx1/2, F_x1, 'g')
texts((Sx1_start + Sx1/2)' + F_x1*1.3 +[0.005 0 0.01]', 'F_{\alpha_1}', 12)

arrow3(Sy1_start + Sy1/2, F_y1, 'g')
texts((Sy1_start + Sy1/2)' + F_y1*1.3 +[0.005 0 0.01]', 'F_{\beta_1}', 12)

arrow3(Sx2_start + Sx2/2, F_x2, 'g')
texts((Sx2_start + Sx2/2)' + F_x2*1.3 +[0.005 0 -0.01]', 'F_{\alpha_2}', 12)

arrow3(Sy2_start + Sy2/2, F_y2, 'g')
texts((Sy2_start + Sy2/2)' + F_y2*1.3 +[0.005 0 -0.01]', 'F_{\beta_2}', 12)
end

tau_y = cross(Sx2,F_y1)
tau_x = cross(Sy1,F_x1)


nice3d
grid off
box off 
axis off

view(30,15);
if(forces)
    zoom(1.5)
end

if(printEPS)
%pause
print('-depsc','-r300', strcat(outputDir,outputName,'.eps'))
end