
if exist('RunFile','var')==0
    clc
    close all
    clear
    WCDTorque
end

%% Init
if exist('RunFile','var')==1

Ngg2 = zeros(1,11^3);
att_w = [0 1 1]';
ii=0;
for A=0:1:10
    %clc; A
    for B=0:1:10
        for C=0:1:10
att_w = [A,B,C]';
att = Crot'*att_w;
att = att/norm(att);
%att = [0.340284153018647 0.647318296998995 0.682045245987906]';
%% Gravity Torque

Ngg=((3*mu)/(ERadiusMean+AltitudeSat)^3)*cross(att,Isat*att);

%disp('Max Gravity Torque')
%norm(Ngg)
ii=ii+1;
Ngg2(ii)=norm(Ngg);
        end
    end
end
Ngg=max(Ngg2)
end