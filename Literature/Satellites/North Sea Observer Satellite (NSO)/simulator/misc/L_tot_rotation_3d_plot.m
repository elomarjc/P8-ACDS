%B field plots'
clc
format long
format compact

size = 12;
X = [1 0 0];
Y = [0 1 0];
Z = [0 0 1];
nul = [0 0 0];

sun=[0.7852904564 0.0761112902 0.6144314206];
sun = [.8 -.9 .88];
sun = [-1 -2 -5]
sun=(1/norm(sun))*sun

e = cross(X,sun);
theta = acos(dot(X,sun)) % Rotation angle
Ltot = [1 1 1]*0.7;

% Calculate Ltotend - rotate ltot theta around e
Ltotend = RotVecArAxe(Ltot,e,-theta);

% Calculate the centre of the circle by projecting Ltot onto e
cc = dot(Ltot,e)/dot(e,e)*e

% Find the intersection points between the circle and R+

% xn = 0 = x*cos(a) + (1 - cos(a))*(c1*c1*x + c1*c2*y + c1*c3*z) + (c2*z - c3*y)*sin(a)
% yn = 0 = y*cos(a) + (1 - cos(a))*(c2*c1*x + c2*c2*y + c2*c3*x) + (c3*x - c1*z)*sin(a)
% zn = 0 = z*cos(a) + (1 - cos(a))*(c3*c1*x + c3*c2*y + c3*c3*z) + (c1*y - c2*x)*sin(a)
% B=A*cphi+(A'*L0)*(1-cphi)*L0+cross(L0,A)*sin(Phi);

en = e/norm(e);
c1 = en(1);
c2 = en(2); 
c3 = en(3); 
x = Ltot(1);
y = Ltot(2);
z = Ltot(3);

%x=0
x11 = (-2*x*c1^2*c2*y*c3*z+c1*c3*z^3*c2^2+c1*c3^3*z*y^2+c1*c2^3*y*z^2+c1*c2*y^3*c3^2-2*c1*c3^2*z^2*y*c2-2*c1*c2^2*y^2*c3*z+c1^2*x*c3^2*y^2+c1^2*x*c2^2*z^2-c1*c2*y*(-(c3*y-c2*z)^2*(-c3^2*y^2+2*c3*y*c2*z+2*x*c1*c2*y+2*x*c1*c3*z-x^2-c2^2*z^2+2*x^2*c1^2))^(1/2)-c1*c3*z*(-(c3*y-c2*z)^2*(-c3^2*y^2+2*c3*y*c2*z+2*x*c1*c2*y+2*x*c1*c3*z-x^2-c2^2*z^2+2*x^2*c1^2))^(1/2)-c1^2*x*(-(c3*y-c2*z)^2*(-c3^2*y^2+2*c3*y*c2*z+2*x*c1*c2*y+2*x*c1*c3*z-x^2-c2^2*z^2+2*x^2*c1^2))^(1/2)+x*(-(c3*y-c2*z)^2*(-c3^2*y^2+2*c3*y*c2*z+2*x*c1*c2*y+2*x*c1*c3*z-x^2-c2^2*z^2+2*x^2*c1^2))^(1/2))/(c3^2*y^2-2*x^2*c1^2+c1^4*x^2+c2^2*z^2-2*x*c1*c2*y-2*x*c1*c3*z+x^2-2*c3*y*c2*z+c1^2*c2^2*y^2+c1^2*c3^2*z^2+2*c1^3*x*c2*y+2*c1^3*x*c3*z+2*c1^2*c2*y*c3*z)/(c3*y-c2*z);
x12 = (-x*c1*c3*z-x*c1*c2*y+c1^2*c3^2*z^2+2*c1^3*x*c2*y+c1^2*c2^2*y^2-x^2*c1^2+c1^4*x^2+2*c1^3*x*c3*z+2*c1^2*c2*y*c3*z+(-(c3*y-c2*z)^2*(-c3^2*y^2+2*c3*y*c2*z+2*x*c1*c2*y+2*x*c1*c3*z-x^2-c2^2*z^2+2*x^2*c1^2))^(1/2))/(c3^2*y^2-2*x^2*c1^2+c1^4*x^2+c2^2*z^2-2*x*c1*c2*y-2*x*c1*c3*z+x^2-2*c3*y*c2*z+c1^2*c2^2*y^2+c1^2*c3^2*z^2+2*c1^3*x*c2*y+2*c1^3*x*c3*z+2*c1^2*c2*y*c3*z);
x21 = (-2*x*c1^2*c2*y*c3*z+c1*c3*z^3*c2^2+c1*c3^3*z*y^2+c1*c2^3*y*z^2+c1*c2*y^3*c3^2-2*c1*c3^2*z^2*y*c2-2*c1*c2^2*y^2*c3*z+c1^2*x*c3^2*y^2+c1^2*x*c2^2*z^2+c1*c2*y*(-(c3*y-c2*z)^2*(-c3^2*y^2+2*c3*y*c2*z+2*x*c1*c2*y+2*x*c1*c3*z-x^2-c2^2*z^2+2*x^2*c1^2))^(1/2)+c1*c3*z*(-(c3*y-c2*z)^2*(-c3^2*y^2+2*c3*y*c2*z+2*x*c1*c2*y+2*x*c1*c3*z-x^2-c2^2*z^2+2*x^2*c1^2))^(1/2)+c1^2*x*(-(c3*y-c2*z)^2*(-c3^2*y^2+2*c3*y*c2*z+2*x*c1*c2*y+2*x*c1*c3*z-x^2-c2^2*z^2+2*x^2*c1^2))^(1/2)-x*(-(c3*y-c2*z)^2*(-c3^2*y^2+2*c3*y*c2*z+2*x*c1*c2*y+2*x*c1*c3*z-x^2-c2^2*z^2+2*x^2*c1^2))^(1/2))/(c3^2*y^2-2*x^2*c1^2+c1^4*x^2+c2^2*z^2-2*x*c1*c2*y-2*x*c1*c3*z+x^2-2*c3*y*c2*z+c1^2*c2^2*y^2+c1^2*c3^2*z^2+2*c1^3*x*c2*y+2*c1^3*x*c3*z+2*c1^2*c2*y*c3*z)/(c3*y-c2*z);
x22 = (-x*c1*c3*z-x*c1*c2*y+c1^2*c3^2*z^2+2*c1^3*x*c2*y+c1^2*c2^2*y^2-x^2*c1^2+c1^4*x^2+2*c1^3*x*c3*z+2*c1^2*c2*y*c3*z-(-(c3*y-c2*z)^2*(-c3^2*y^2+2*c3*y*c2*z+2*x*c1*c2*y+2*x*c1*c3*z-x^2-c2^2*z^2+2*x^2*c1^2))^(1/2))/(c3^2*y^2-2*x^2*c1^2+c1^4*x^2+c2^2*z^2-2*x*c1*c2*y-2*x*c1*c3*z+x^2-2*c3*y*c2*z+c1^2*c2^2*y^2+c1^2*c3^2*z^2+2*c1^3*x*c2*y+2*c1^3*x*c3*z+2*c1^2*c2*y*c3*z);
if(isreal([x11 x12 x21 x22]))
    alphax1 = atan2(x11,x12)
    alphax2 = atan2(x21,x22)
    Ltotx1 = RotVecArAxe(Ltot,e,alphax1);
    Ltotx2 = RotVecArAxe(Ltot,e,alphax2);
else
    Ltotx1 = nul';
    Ltotx2 = nul';
end


%y=0
y11 = -(c2*c3^3*x^3+x*c1^3*c2*z^2-2*c2*c3^2*x^2*c1*z+x^3*c1*c2*c3^2-c2*c3*x*((c3*x-c1*z)^2*(x^2*c3^2-2*y*c2*c3*x-2*x*c1*c3*z-2*x*c1*c2*y+y^2-2*y^2*c2^2+c1^2*z^2))^(1/2)+c2*c3*x*c1^2*z^2-2*x^2*c1^2*c3*c2*z-2*c2^2*y*x*c1*c3*z+c2^2*y*c1^2*z^2+y*c3^2*c2^2*x^2-c2^2*y*((c3*x-c1*z)^2*(x^2*c3^2-2*y*c2*c3*x-2*x*c1*c3*z-2*x*c1*c2*y+y^2-2*y^2*c2^2+c1^2*z^2))^(1/2)-x*c1*c2*((c3*x-c1*z)^2*(x^2*c3^2-2*y*c2*c3*x-2*x*c1*c3*z-2*x*c1*c2*y+y^2-2*y^2*c2^2+c1^2*z^2))^(1/2)+y*((c3*x-c1*z)^2*(x^2*c3^2-2*y*c2*c3*x-2*x*c1*c3*z-2*x*c1*c2*y+y^2-2*y^2*c2^2+c1^2*z^2))^(1/2))/(-2*x*c1*c3*z-2*x*c1*c2*y-2*y*c2*c3*x+x^2*c1^2*c2^2+c3^2*c2^2*x^2-2*y^2*c2^2+2*x*c1*c2^3*y+2*x^2*c1*c2^2*c3+c2^4*y^2+x^2*c3^2+c1^2*z^2+2*c2^3*y*c3*x+y^2)/(c3*x-c1*z);
y12 = (-x*c1*c2*y+2*x*c1*c2^3*y+2*c2^3*y*c3*x-y^2*c2^2+x^2*c1^2*c2^2+c3^2*c2^2*x^2+2*x^2*c1*c2^2*c3-y*c2*c3*x+c2^4*y^2+((c3*x-c1*z)^2*(x^2*c3^2-2*y*c2*c3*x-2*x*c1*c3*z-2*x*c1*c2*y+y^2-2*y^2*c2^2+c1^2*z^2))^(1/2))/(-2*x*c1*c3*z-2*x*c1*c2*y-2*y*c2*c3*x+x^2*c1^2*c2^2+c3^2*c2^2*x^2-2*y^2*c2^2+2*x*c1*c2^3*y+2*x^2*c1*c2^2*c3+c2^4*y^2+x^2*c3^2+c1^2*z^2+2*c2^3*y*c3*x+y^2);
y21 = -(c2*c3^3*x^3+x*c1^3*c2*z^2-2*c2*c3^2*x^2*c1*z+x^3*c1*c2*c3^2+c2*c3*x*((c3*x-c1*z)^2*(x^2*c3^2-2*y*c2*c3*x-2*x*c1*c3*z-2*x*c1*c2*y+y^2-2*y^2*c2^2+c1^2*z^2))^(1/2)+c2*c3*x*c1^2*z^2-2*x^2*c1^2*c3*c2*z-2*c2^2*y*x*c1*c3*z+c2^2*y*c1^2*z^2+y*c3^2*c2^2*x^2+c2^2*y*((c3*x-c1*z)^2*(x^2*c3^2-2*y*c2*c3*x-2*x*c1*c3*z-2*x*c1*c2*y+y^2-2*y^2*c2^2+c1^2*z^2))^(1/2)+x*c1*c2*((c3*x-c1*z)^2*(x^2*c3^2-2*y*c2*c3*x-2*x*c1*c3*z-2*x*c1*c2*y+y^2-2*y^2*c2^2+c1^2*z^2))^(1/2)-y*((c3*x-c1*z)^2*(x^2*c3^2-2*y*c2*c3*x-2*x*c1*c3*z-2*x*c1*c2*y+y^2-2*y^2*c2^2+c1^2*z^2))^(1/2))/(-2*x*c1*c3*z-2*x*c1*c2*y-2*y*c2*c3*x+x^2*c1^2*c2^2+c3^2*c2^2*x^2-2*y^2*c2^2+2*x*c1*c2^3*y+2*x^2*c1*c2^2*c3+c2^4*y^2+x^2*c3^2+c1^2*z^2+2*c2^3*y*c3*x+y^2)/(c3*x-c1*z);
y22 = (-x*c1*c2*y+2*x*c1*c2^3*y+2*c2^3*y*c3*x-y^2*c2^2+x^2*c1^2*c2^2+c3^2*c2^2*x^2+2*x^2*c1*c2^2*c3-y*c2*c3*x+c2^4*y^2-((c3*x-c1*z)^2*(x^2*c3^2-2*y*c2*c3*x-2*x*c1*c3*z-2*x*c1*c2*y+y^2-2*y^2*c2^2+c1^2*z^2))^(1/2))/(-2*x*c1*c3*z-2*x*c1*c2*y-2*y*c2*c3*x+x^2*c1^2*c2^2+c3^2*c2^2*x^2-2*y^2*c2^2+2*x*c1*c2^3*y+2*x^2*c1*c2^2*c3+c2^4*y^2+x^2*c3^2+c1^2*z^2+2*c2^3*y*c3*x+y^2);
if(isreal([y11 y12 y21 y22]))
    alphay1 = atan2(y11,y12)
    alphay2 = atan2(y21,y22)
    Ltoty1 = RotVecArAxe(Ltot,e,alphay1);
    Ltoty2 = RotVecArAxe(Ltot,e,alphay2);
else
    Ltoty1 = nul';
    Ltoty2 = nul';
end

%z=0
z11 = (c3^2*z*c2^2*x^2+c3^2*z*c1^2*y^2-c3*c1*x*((c2*x-c1*y)^2*(c2^2*x^2-2*x*c1*c2*y-2*x*c1*c3*z-2*c3*y*c2*z-2*z^2*c3^2+z^2+c1^2*y^2))^(1/2)-c3*c2*y*((c2*x-c1*y)^2*(c2^2*x^2-2*x*c1*c2*y-2*x*c1*c3*z-2*c3*y*c2*z-2*z^2*c3^2+z^2+c1^2*y^2))^(1/2)-c3^2*z*((c2*x-c1*y)^2*(c2^2*x^2-2*x*c1*c2*y-2*x*c1*c3*z-2*c3*y*c2*z-2*z^2*c3^2+z^2+c1^2*y^2))^(1/2)-2*x*c1*c2^2*y^2*c3-2*x^2*c1^2*c3*y*c2+c3*y*c2^3*x^2+c3*y^3*c2*c1^2+x*c1^3*c3*y^2+x^3*c1*c3*c2^2+z*((c2*x-c1*y)^2*(c2^2*x^2-2*x*c1*c2*y-2*x*c1*c3*z-2*c3*y*c2*z-2*z^2*c3^2+z^2+c1^2*y^2))^(1/2)-2*z*c3^2*c1*x*c2*y)/(-2*z^2*c3^2+c3^4*z^2+2*c3^2*c1*x*c2*y+2*x*c1*c3^3*z+2*c3^3*c2*y*z-2*x*c1*c2*y+c2^2*x^2+c1^2*y^2-2*x*c1*c3*z-2*c3*y*c2*z+x^2*c1^2*c3^2+c3^2*c2^2*y^2+z^2)/(c2*x-c1*y);
z12 = (-z^2*c3^2-c3*y*c2*z+2*x*c1*c3^3*z+c3^4*z^2-x*c1*c3*z+2*c3^2*c1*x*c2*y+c3^2*c2^2*y^2+2*c3^3*c2*y*z+x^2*c1^2*c3^2+((c2*x-c1*y)^2*(c2^2*x^2-2*x*c1*c2*y-2*x*c1*c3*z-2*c3*y*c2*z-2*z^2*c3^2+z^2+c1^2*y^2))^(1/2))/(-2*z^2*c3^2+c3^4*z^2+2*c3^2*c1*x*c2*y+2*x*c1*c3^3*z+2*c3^3*c2*y*z-2*x*c1*c2*y+c2^2*x^2+c1^2*y^2-2*x*c1*c3*z-2*c3*y*c2*z+x^2*c1^2*c3^2+c3^2*c2^2*y^2+z^2);
z21 = (c3^2*z*c2^2*x^2+c3^2*z*c1^2*y^2+c3*c1*x*((c2*x-c1*y)^2*(c2^2*x^2-2*x*c1*c2*y-2*x*c1*c3*z-2*c3*y*c2*z-2*z^2*c3^2+z^2+c1^2*y^2))^(1/2)+c3*c2*y*((c2*x-c1*y)^2*(c2^2*x^2-2*x*c1*c2*y-2*x*c1*c3*z-2*c3*y*c2*z-2*z^2*c3^2+z^2+c1^2*y^2))^(1/2)+c3^2*z*((c2*x-c1*y)^2*(c2^2*x^2-2*x*c1*c2*y-2*x*c1*c3*z-2*c3*y*c2*z-2*z^2*c3^2+z^2+c1^2*y^2))^(1/2)-2*x*c1*c2^2*y^2*c3-2*x^2*c1^2*c3*y*c2+c3*y*c2^3*x^2+c3*y^3*c2*c1^2+x*c1^3*c3*y^2+x^3*c1*c3*c2^2-z*((c2*x-c1*y)^2*(c2^2*x^2-2*x*c1*c2*y-2*x*c1*c3*z-2*c3*y*c2*z-2*z^2*c3^2+z^2+c1^2*y^2))^(1/2)-2*z*c3^2*c1*x*c2*y)/(-2*z^2*c3^2+c3^4*z^2+2*c3^2*c1*x*c2*y+2*x*c1*c3^3*z+2*c3^3*c2*y*z-2*x*c1*c2*y+c2^2*x^2+c1^2*y^2-2*x*c1*c3*z-2*c3*y*c2*z+x^2*c1^2*c3^2+c3^2*c2^2*y^2+z^2)/(c2*x-c1*y);
z22 = (-z^2*c3^2-c3*y*c2*z+2*x*c1*c3^3*z+c3^4*z^2-x*c1*c3*z+2*c3^2*c1*x*c2*y+c3^2*c2^2*y^2+2*c3^3*c2*y*z+x^2*c1^2*c3^2-((c2*x-c1*y)^2*(c2^2*x^2-2*x*c1*c2*y-2*x*c1*c3*z-2*c3*y*c2*z-2*z^2*c3^2+z^2+c1^2*y^2))^(1/2))/(-2*z^2*c3^2+c3^4*z^2+2*c3^2*c1*x*c2*y+2*x*c1*c3^3*z+2*c3^3*c2*y*z-2*x*c1*c2*y+c2^2*x^2+c1^2*y^2-2*x*c1*c3*z-2*c3*y*c2*z+x^2*c1^2*c3^2+c3^2*c2^2*y^2+z^2);
if(isreal([z11 z12 z21 z22]))
    alphaz1 = atan2(z11,z12)
    alphaz2 = atan2(z21,z22)
    Ltotz1 = RotVecArAxe(Ltot,e,alphaz1);
    Ltotz2 = RotVecArAxe(Ltot,e,alphaz2);
else
    Ltotz1 = nul';
    Ltotz2 = nul';
end

% Calculate the circle plot points
theta_c = -(0:.1:2*pi);
Ltotcirc = [0 0 0]';
for n = 1:length(theta_c)
    Ltotcirc(:,n) = RotVecArAxe(Ltot,e,theta_c(n));
end

% Calculate Xend - rotate X theta around e
Xend = RotVecArAxe([1 0 0],e,theta);

% Calculate the circle plot points
theta_c = 0:.1:2*pi;
Xcirc = [0 0 0]';
for n = 1:length(theta_c)
    Xcirc(:,n) = RotVecArAxe([1 0 0],e,theta_c(n));
end

dir = -1
theta = theta*dir

intangles = [];
if(dir == sign(alphax1))
    intangles = [intangles alphax1]
end
if(dir == sign(alphax2))
    intangles = [intangles alphax2]
end
if(dir == sign(alphay1))
    intangles = [intangles alphay1]
end
if(dir == sign(alphay2))
    intangles = [intangles alphay2]
end
if(dir == sign(alphaz1))
    intangles = [intangles alphaz1]
end
if(dir == sign(alphaz2))
    intangles = [intangles alphaz2]
end

intangle = dir*min(abs(intangles))

intvect = RotVecArAxe(Ltot,e,intangle);


% % Determine which intersection vector to use
% % Find the direction of rotation (round of for calculation errors)
% dir = sign(round(cross(Ltot,Ltotend)*1E6)/1E6)
% ix1 = sign(round(cross(Ltot,Ltotx1)*1E6)/1E6)
% ix2 = sign(round(cross(Ltot,Ltotx2)*1E6)/1E6)
% iy1 = sign(round(cross(Ltot,Ltoty1)*1E6)/1E6)
% iy2 = sign(round(cross(Ltot,Ltoty2)*1E6)/1E6)
% iz1 = sign(round(cross(Ltot,Ltotz1)*1E6)/1E6)
% iz2 = sign(round(cross(Ltot,Ltotz2)*1E6)/1E6)
% 
% % % 
% ixe1 = sign(round(cross(Ltotend,ix1)*1E6)/1E6)
% ixe2 = sign(round(cross(Ltotend,ix2)*1E6)/1E6)
% iye1 = sign(round(cross(Ltotend,iy1)*1E6)/1E6)
% iye2 = sign(round(cross(Ltotend,iy2)*1E6)/1E6)
% ize1 = sign(round(cross(Ltotend,iz1)*1E6)/1E6)
% ize2 = sign(round(cross(Ltotend,iz2)*1E6)/1E6)
% 
% % Find the intersection vector by checking that the direction is right and
% % that the vector lies between Ltotp and Ltotendp
% intvect = nul;
% 
% if((dir == ix1) & (ix1 == -ixe1))
%     intvect = Ltotx1
% end
% if((dir == ix2) & (ix2 == -ixe2))
%     intvect = Ltotx2
% end
% if((dir == iy1) & (iy1 == -iye1))
%     intvect = Ltoty1
% end
% if((dir == iz2) & (iz2 == -ize2))
%     intvect = Ltotz2
% end
% if((dir == iz1) & (iz1 == -ize1))
%     intvect = Ltotz1
% end
% if((dir == iz2) & (iz2 == -ize2))
%     intvect = Ltotz2
% end
% Normalize the intersection vector
%intvect = intvect/norm(intvect);

% Find the maximum alowable turning angle
% Check if Ltotend is in R+




figure(1)
clf
arrow3(nul, [1 0 0], 'k'); nice3d
hold on
texts([1.1 0 0], 'x',size)
arrow3(nul, [0 1 0], 'k'); nice3d
texts([0 1.1 0], 'y',size)
arrow3(nul, [0 0 1], 'k'); nice3d
texts([0 0 1.1], 'z',size)

% Plot the e vector
arrow3(nul, e, 'b'); nice3d
texts(e*1.1, 'e',size)

% Plot the Ltot vector
arrow3(nul, Ltot, 'c'); nice3d
texts(Ltot*1.1, 'Ltot',size)

% Plot cc
arrow3(nul, cc, 'y'); nice3d
texts(cc*1.1, 'CC',size)

% Plot Ltot intersection vectors
% X
if(max(abs(Ltotx1))>0)
    disp('janx');
    arrow3(nul, Ltotx1, 'g'); nice3d
    texts(Ltotx1*1.1, 'Ltotx1',size)
    arrow3(nul, Ltotx2, 'g'); nice3d
    texts(Ltotx2*1.1, 'Ltotx2',size)
end

% Y
if(max(abs(Ltotx1))>0')
    disp('jany');
    arrow3(nul, Ltoty1, 'g'); nice3d
    texts(Ltoty1*1.1, 'Ltoty1',size)
    arrow3(nul, Ltoty2, 'g'); nice3d
    texts(Ltoty2*1.1, 'Ltoty2',size)
end

% Z
if(max(abs(Ltotx1))>0')
    disp('janz');
    arrow3(nul, Ltotz1, 'g'); nice3d
    texts(Ltotz1*1.1, 'Ltotz1',size)
    arrow3(nul, Ltotz2, 'g'); nice3d
    texts(Ltotz2*1.1, 'Ltotz2',size)
end

% Plot Ltotend
arrow3(nul, Ltotend, 'c'); nice3d
texts(Ltotend*1.1, 'Ltotend',size)

% % Plot intersection vector
arrow3(nul, intvect*1.5, 'r'); nice3d
texts(intvect*1.6, 'intvect',size)

% Plot the circle
for n = 2:length(Ltotcirc)
    x = [Ltotcirc(1,n-1) Ltotcirc(1,n)];
    y = [Ltotcirc(2,n-1) Ltotcirc(2,n)];
    z = [Ltotcirc(3,n-1) Ltotcirc(3,n)];
    plot3(x,y,z,'c','linewidth',2)
end

arrow3(nul, sun, 'r'); nice3d
texts(sun*1.1, 'Sun',size)

% Plot the circle
for n = 2:length(Xcirc)
    x = [Xcirc(1,n-1) Xcirc(1,n)];
    y = [Xcirc(2,n-1) Xcirc(2,n)];
    z = [Xcirc(3,n-1) Xcirc(3,n)];
    plot3(x,y,z,'r','linewidth',2)
end


nice3d;
axis off
grid on