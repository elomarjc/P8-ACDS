%% Clear, clean and close all
clear
clc
close all

%% Init

% [m]
SatLength=0.2;
SatHeigth=0.1;
SatWidth=0.1;

% [m]
AntLength=0.39;
AntRadius=0.08;

% [m]
BoardLength=0.87;
BoardWidth=0.87;

% [kg]
SatWeight=1.64;

% Precission
prec=(pi/2)*0.01;
TotalLength=SatLength+(AntLength*2)

%% A1 (0.1 X 0.2)

A1=zeros(1,length((pi/2)/prec)+1);
A2=A1;
A3=A1;
for i=1:(((pi/2)/prec)+2)
    j=1;
    for j=1:(((pi/2)/prec)+2)
        A1(j)=cos((j-1)*prec)*cos((i-1)*prec)*SatLength*SatWidth;
        A2(j)=sin((j-1)*prec)*cos((i-1)*prec)*SatWidth*SatHeigth;
        A3(j)=sin((i-1)*prec)*cos((j-1)*prec)*SatHeigth*SatLength;
        Atot1(j)=A1(j)+A2(j)+A3(j);
    end
    A11(i,1:j)=A1;
    A22(i,1:j)=A2;
    A33(i,1:j)=A3;
    Atot(i,1:j)=Atot1;
    Asummax(i)=max(Atot1);
end

Amax=max(Asummax)



