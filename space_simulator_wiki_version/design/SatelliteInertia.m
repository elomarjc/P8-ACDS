if exist('RunFile','var')==0
    clc
    close all
    clear
    WCDTorque
end

%% Init
if exist('RunFile','var')==1

Extra_batt=0; % 0=2 batt's, 1=4 batt's

SideDepth=0.5e-3;

BoardLength=87e-3;
BoardWidth=87e-3;
BoardDepth=2e-3;
GPSLength=70e-3;
GPSWidth=46e-3;
AntennaWidth=100e-3;
AntennaLength=100e-3;
BattLength=65e-3;
BattDia=15e-3;

AntennaHeigth=20e-3;
GPSHeigth=12e-3;
ADCSHeigth=5e-3;
FPGAHeigth=14e-3;
AISHeigth=12e-3;
UHFHeigth=12e-3;
EPSHeigth=12e-3;
BattHeigth=8e-3;

mPlate=10.67e-3; % Weigth of one side plate
mPCB=3e-3; % Weigth of small PCB for sun sensor
mCoil=17.76e-3; % Weigth of one set of magnetorquers

mA1=(mPlate+mPCB+mCoil); %Plate (10x10x0.05) + torquer + small PCB
mB1=mA1;
mB2=mA1;
mA2=(mPlate+mPCB); %Plate (10x10x0.05) + small PCB
mB3=mA2;
mB4=mA2;

factor_nylon2alu=2.6; % Factor between the density of nylon and aluminium.
mAntenna=77.6*2*factor_nylon2alu*1e-3; % 2 Antenna modules in teflon
mGPS=22.48e-3;
mADCS=50e-3;
mFPGA=30e-3;    % Guestimated
mAIS=(32.17+51.17-23.57)*1e-3;    % Combined weigth of old AIS1 & 2
mUHF=(33.18+4)*1e-3;
mBatt=41.57e-3;  % One batt
mEPS=43.16e-3; % EPS

% CoM displacement of 2 cm from GoM (3.2 cm from guess):
% 4xBatt, Batt=155e-3, Antenna=77.6*2*1e-3

% CoM displacement of -2 cm from GoM (0.8 cm from guess):
% 2xBatt, Antenna=77.6*2*4.6*1e-3

FrameADepth=4.5e-3;
FrameBDepth=8.5e-3;
mFramePart=(93.27/12)*1e-3;

% Vector with heights.
BoardHeigthVec=[AntennaHeigth,GPSHeigth,ADCSHeigth,FPGAHeigth,AISHeigth,UHFHeigth,EPSHeigth,BattDia];

% Vector with weights.
SideWeigthVec=[mA1,mA2,mB1,mB2,mB3,mB4];
BoardWeigthVec=[mAntenna,mGPS,mADCS,mFPGA,mAIS,mUHF,mEPS,mBatt,mBatt];
if Extra_batt==1
    BoardWeigthVec=[BoardWeigthVec,mBatt,mBatt];
end
FrameWeigthVec=ones(1,12)*mFramePart;
mVec=[SideWeigthVec,FrameWeigthVec,BoardWeigthVec];
mSat=sum(mVec);

HeigthTotal=sum(BoardHeigthVec);
%--------    0 mm
% Antenna    10 mm
%--------    20 mm
% GPS       
%--------    25 mm
% ADCS       37 mm
%--------
% FPGA       51 mm
%--------
% AIS        63 mm
%--------
% UHF        75 mm
%--------
% EPS        87 mm
% Batt 
%--------   110 mm
%% Inertia of elements

% Plates
IsideA1=[(1/12)*mA1*(SideDepth^2+SatHeigth^2) 0 0;
         0 (1/12)*mA1*(SatWidth^2+SatHeigth^2) 0;
         0 0 (1/12)*mA1*(SideDepth^2+SatWidth^2)];

IsideA2=[(1/12)*mA2*(SideDepth^2+SatHeigth^2) 0 0;
         0 (1/12)*mA2*(SatWidth^2+SatHeigth^2) 0;
         0 0 (1/12)*mA2*(SideDepth^2+SatWidth^2)];

IsideB1=[(1/12)*mB1*(SatHeigth^2+SatLength^2) 0 0;
         0 (1/12)*mB1*(SatHeigth^2+SideDepth^2) 0;
         0 0 (1/12)*mB1*(SideDepth^2+SatLength^2)];

IsideB3=[(1/12)*mB3*(SatHeigth^2+SatLength^2) 0 0;
         0 (1/12)*mB3*(SatHeigth^2+SideDepth^2) 0;
         0 0 (1/12)*mB3*(SideDepth^2+SatLength^2)];

IsideB2=[(1/12)*mB2*(SideDepth^2+SatLength^2) 0 0;
         0 (1/12)*mB2*(SatWidth^2+SideDepth^2) 0;
         0 0 (1/12)*mB2*(SatWidth^2+SatLength^2)];

IsideB4=[(1/12)*mB4*(SideDepth^2+SatLength^2) 0 0;
         0 (1/12)*mB4*(SatWidth^2+SideDepth^2) 0;
         0 0 (1/12)*mB4*(SatWidth^2+SatLength^2)];

% Frame

IframeA1=[(1/12)*mFramePart*(FrameBDepth^2+SatHeigth^2) 0 0;
         0 (1/12)*mFramePart*(FrameBDepth^2+SatHeigth^2) 0;
         0 0 (1/12)*mFramePart*(FrameBDepth^2+FrameBDepth^2)];

IframeA3=IframeA1;
IframeA5=IframeA1;
IframeA7=IframeA1;

IframeA4=[(1/12)*mFramePart*(FrameADepth^2+FrameADepth^2) 0 0;
         0 (1/12)*mFramePart*(FrameADepth^2+SatWidth^2) 0;
         0 0 (1/12)*mFramePart*(SatWidth^2+FrameADepth^2)];

IframeA2=IframeA4;
IframeA6=IframeA4;
IframeA8=IframeA4;

IframeB1=[(1/12)*mFramePart*(FrameADepth^2+SatLength^2) 0 0;
         0 (1/12)*mFramePart*(FrameADepth^2+FrameADepth^2) 0;
         0 0 (1/12)*mFramePart*(FrameADepth^2+SatLength^2)];

IframeB2=IframeB1;
IframeB3=IframeB1;
IframeB4=IframeB1;

% Boards

IAntenna=[(1/12)*mAntenna*(AntennaLength^2+AntennaHeigth^2) 0 0;
       0 (1/12)*mAntenna*(AntennaHeigth^2+AntennaWidth^2) 0;
       0 0 (1/12)*mAntenna*(AntennaWidth^2+AntennaLength^2)];

IGPS=[(1/12)*mGPS*(GPSLength^2+GPSHeigth^2) 0 0;
       0 (1/12)*mGPS*(GPSHeigth^2+GPSWidth^2) 0;
       0 0 (1/12)*mGPS*(BoardWidth^2+GPSLength^2)];

IADCS=[(1/12)*mADCS*(BoardLength^2+BoardDepth^2) 0 0;
       0 (1/12)*mADCS*(BoardDepth^2+BoardWidth^2) 0;
       0 0 (1/12)*mADCS*(BoardWidth^2+BoardLength^2)];

IFPGA=[(1/12)*mFPGA*(BoardLength^2+BoardDepth^2) 0 0;
       0 (1/12)*mFPGA*(BoardDepth^2+BoardWidth^2) 0;
       0 0 (1/12)*mFPGA*(BoardWidth^2+BoardLength^2)];
   
IAIS=[(1/12)*mAIS*(BoardLength^2+BoardDepth^2) 0 0;
       0 (1/12)*mAIS*(BoardDepth^2+BoardWidth^2) 0;
       0 0 (1/12)*mAIS*(BoardWidth^2+BoardLength^2)];

IUHF=[(1/12)*mUHF*(BoardLength^2+BoardDepth^2) 0 0;
       0 (1/12)*mUHF*(BoardDepth^2+BoardWidth^2) 0;
       0 0 (1/12)*mUHF*(BoardWidth^2+BoardLength^2)];

IEPS=[(1/12)*mEPS*(BoardLength^2+EPSHeigth^2) 0 0;
       0 (1/12)*mEPS*(EPSHeigth^2+BoardWidth^2) 0;
       0 0 (1/12)*mEPS*(BoardWidth^2+BoardLength^2)];

IBatt1=[(1/12)*mBatt*(3*(BattDia/2)^2+BattLength^2) 0 0;
       0 (1/2)*mBatt*(BattDia/2)^2 0;
       0 0 (1/12)*mBatt*(3*(BattDia/2)^2+BattLength^2)];

IBatt2=IBatt1;
IBatt3=IBatt1;
IBatt4=IBatt1;

%% CoM Sat

% Plates
pA1=[SatWidth/2,0,SatHeigth/2]';
pA2=[SatWidth/2,SatLength,SatHeigth/2]';
pB1=[0,SatLength/2,SatHeigth/2]';
pB2=[SatWidth/2,SatLength/2,SatHeigth]';
pB3=[SatWidth,SatLength/2,SatHeigth/2]';
pB4=[SatWidth/2,SatLength/2,0]';

% Frame
pFA1=[0,0,SatHeigth/2]';
pFA2=[SatWidth/2,0,AntennaHeigth+(FrameADepth/2)]';
pFA3=[SatWidth,0,SatHeigth/2]';
pFA4=[SatWidth/2,0,SatHeigth]';
pFA5=[0,SatLength,SatHeigth/2]';
pFA6=[SatWidth/2,SatLength,AntennaHeigth+(FrameADepth/2)]';
pFA7=[SatWidth,SatLength,SatHeigth/2]';
pFA8=[SatWidth/2,SatLength,SatHeigth]';
pFB1=[0,SatLength/2,SatHeigth]';
pFB2=[0,SatLength/2,AntennaHeigth+(FrameADepth/2)]';
pFB3=[SatWidth,SatLength/2,AntennaHeigth+(FrameADepth/2)]';
pFB4=[SatWidth,SatLength/2,SatHeigth]';

% Boards
pAntenna=[SatWidth/2,SatLength/2,AntennaHeigth/2]';
pGPS=[SatWidth/2,(SatLength/2)-7e-3,AntennaHeigth+GPSHeigth]';
pADCS=[SatWidth/2,SatLength/2,AntennaHeigth+GPSHeigth+ADCSHeigth]';
pFPGA=[SatWidth/2,SatLength/2,AntennaHeigth+GPSHeigth+ADCSHeigth+FPGAHeigth]';
pAIS=[SatWidth/2,SatLength/2,AntennaHeigth+GPSHeigth+ADCSHeigth+FPGAHeigth+AISHeigth]';
pUHF=[SatWidth/2,SatLength/2,AntennaHeigth+GPSHeigth+ADCSHeigth+FPGAHeigth+AISHeigth+UHFHeigth]';
pEPS=[SatWidth/2,SatLength/2,AntennaHeigth+GPSHeigth+ADCSHeigth+FPGAHeigth+AISHeigth+UHFHeigth+EPSHeigth]';
pBatt1=[(SatWidth/2)-BattDia/2,SatLength/2,AntennaHeigth+GPSHeigth+ADCSHeigth+FPGAHeigth+AISHeigth+UHFHeigth+EPSHeigth+BattDia]';
pBatt2=[(SatWidth/2)+BattDia/2,SatLength/2,AntennaHeigth+GPSHeigth+ADCSHeigth+FPGAHeigth+AISHeigth+UHFHeigth+EPSHeigth+BattDia]';
pBatt3=[(SatWidth/2)-(BattDia/2)*3,SatLength/2,AntennaHeigth+GPSHeigth+ADCSHeigth+FPGAHeigth+AISHeigth+UHFHeigth+EPSHeigth+BattDia]';
pBatt4=[(SatWidth/2)+(BattDia/2)*3,SatLength/2,AntennaHeigth+GPSHeigth+ADCSHeigth+FPGAHeigth+AISHeigth+UHFHeigth+EPSHeigth+BattDia]';

pVec=[pA1,pA2,pB1,pB2,pB3,pB4,pFA1,pFA2,pFA3,pFA4,pFA5,pFA6,pFA7,pFA8,pFB1,pFB2,pFB3,pFB4,pAntenna,pGPS,pADCS,pFPGA,pAIS,pUHF,pEPS,pBatt1,pBatt2];

if Extra_batt==1
    pVec=[pVec,pBatt3,pBatt4];
end

pCoM=[0;0;0];
for i=1:length(pVec)
    pTemp(:,1)=mVec(i)*pVec(:,i);
    pCoM=pCoM+pTemp;
end

pCoM=pCoM/sum(mVec);

%if MaxTorqueEnable==1
%    pCoM(2)=pCoM(2)-20e-3;
%    CoM=pCoM;
%else
    CoM=pCoM;
%end


%% Inertia Sat

IA1=IsideA1+mA1*(((pCoM-pA1)'*(pCoM-pA1)*eye(3))-((pCoM-pA1)*(pCoM-pA1)'));
IA2=IsideA2+mA2*(((pCoM-pA2)'*(pCoM-pA2)*eye(3))-((pCoM-pA2)*(pCoM-pA2)'));
IB1=IsideB1+mB1*(((pCoM-pB1)'*(pCoM-pB1)*eye(3))-((pCoM-pB1)*(pCoM-pB1)'));
IB2=IsideB2+mB2*(((pCoM-pB2)'*(pCoM-pB2)*eye(3))-((pCoM-pB2)*(pCoM-pB2)'));
IB3=IsideB3+mB3*(((pCoM-pB3)'*(pCoM-pB3)*eye(3))-((pCoM-pB3)*(pCoM-pB3)'));
IB4=IsideB4+mB4*(((pCoM-pB4)'*(pCoM-pB4)*eye(3))-((pCoM-pB4)*(pCoM-pB4)'));

IFA1=IframeA1+mFramePart*(((pCoM-pFA1)'*(pCoM-pFA1)*eye(3))-((pCoM-pFA1)*(pCoM-pFA1)'));
IFA2=IframeA2+mFramePart*(((pCoM-pFA2)'*(pCoM-pFA2)*eye(3))-((pCoM-pFA2)*(pCoM-pFA2)'));
IFA3=IframeA3+mFramePart*(((pCoM-pFA3)'*(pCoM-pFA3)*eye(3))-((pCoM-pFA3)*(pCoM-pFA3)'));
IFA4=IframeA4+mFramePart*(((pCoM-pFA4)'*(pCoM-pFA4)*eye(3))-((pCoM-pFA4)*(pCoM-pFA4)'));
IFA5=IframeA5+mFramePart*(((pCoM-pFA5)'*(pCoM-pFA5)*eye(3))-((pCoM-pFA5)*(pCoM-pFA5)'));
IFA6=IframeA6+mFramePart*(((pCoM-pFA6)'*(pCoM-pFA6)*eye(3))-((pCoM-pFA6)*(pCoM-pFA6)'));
IFA7=IframeA7+mFramePart*(((pCoM-pFA7)'*(pCoM-pFA7)*eye(3))-((pCoM-pFA7)*(pCoM-pFA7)'));
IFA8=IframeA8+mFramePart*(((pCoM-pFA8)'*(pCoM-pFA8)*eye(3))-((pCoM-pFA8)*(pCoM-pFA8)'));
IFB1=IframeB1+mFramePart*(((pCoM-pFB1)'*(pCoM-pFB1)*eye(3))-((pCoM-pFB1)*(pCoM-pFB1)'));
IFB2=IframeB2+mFramePart*(((pCoM-pFB2)'*(pCoM-pFB2)*eye(3))-((pCoM-pFB2)*(pCoM-pFB2)'));
IFB3=IframeB3+mFramePart*(((pCoM-pFB3)'*(pCoM-pFB3)*eye(3))-((pCoM-pFB3)*(pCoM-pFB3)'));
IFB4=IframeB4+mFramePart*(((pCoM-pFB4)'*(pCoM-pFB4)*eye(3))-((pCoM-pFB4)*(pCoM-pFB4)'));

IAntenna=IAntenna+mAntenna*(((pCoM-pAntenna)'*(pCoM-pAntenna)*eye(3))-((pCoM-pAntenna)*(pCoM-pAntenna)'));
IGPS=IGPS+mGPS*(((pCoM-pGPS)'*(pCoM-pGPS)*eye(3))-((pCoM-pGPS)*(pCoM-pGPS)'));
IADCS=IADCS+mADCS*(((pCoM-pADCS)'*(pCoM-pADCS)*eye(3))-((pCoM-pADCS)*(pCoM-pADCS)'));
IFPGA=IFPGA+mFPGA*(((pCoM-pFPGA)'*(pCoM-pFPGA)*eye(3))-((pCoM-pFPGA)*(pCoM-pFPGA)'));
IAIS=IAIS+mAIS*(((pCoM-pAIS)'*(pCoM-pAIS)*eye(3))-((pCoM-pAIS)*(pCoM-pAIS)'));
IUHF=IUHF+mUHF*(((pCoM-pUHF)'*(pCoM-pUHF)*eye(3))-((pCoM-pUHF)*(pCoM-pUHF)'));
IEPS=IEPS+mEPS*(((pCoM-pEPS)'*(pCoM-pEPS)*eye(3))-((pCoM-pEPS)*(pCoM-pEPS)'));
IBatt1=IBatt1+mBatt*(((pCoM-pBatt1)'*(pCoM-pBatt1)*eye(3))-((pCoM-pBatt1)*(pCoM-pBatt1)'));
IBatt2=IBatt2+mBatt*(((pCoM-pBatt2)'*(pCoM-pBatt2)*eye(3))-((pCoM-pBatt2)*(pCoM-pBatt2)'));
IBatt3=IBatt3+mBatt*(((pCoM-pBatt3)'*(pCoM-pBatt3)*eye(3))-((pCoM-pBatt3)*(pCoM-pBatt3)'));
IBatt4=IBatt4+mBatt*(((pCoM-pBatt4)'*(pCoM-pBatt4)*eye(3))-((pCoM-pBatt4)*(pCoM-pBatt4)'));

Isat=(IA1+IA2+IB1+IB2+IB3+IB4)+(IFA1+IFA2+IFA3+IFA4+IFA5+IFA6+IFA7+IFA8+IFB1+IFB2+IFB3+IFB4)+(IAntenna+IGPS+IADCS+IFPGA+IAIS+IUHF+IEPS+IBatt1+IBatt2)

if Extra_batt==1
    Isat=Isat+IBatt3+IBatt4;
end

[A,Isat]=eig(Isat);
Crot=A;

    % Attitude matrix to quaternion conversion
    % http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixT
    % oQuaternion/index.htm
    Trace = A(1,1) + A(2,2) + A(3,3);
    if Trace > 0
        S  = sqrt(Trace+1.0) * 2;
        q4 = 0.25 * S;
        q1 = (A(3,2) - A(2,3)) / S;
        q2 = (A(1,3) - A(3,1)) / S;  
        q3 = (A(2,1) - A(1,2)) / S;
    elseif (A(1,1) > A(2,2)) && (A(1,1) > A(3,3))
        S  = sqrt(1.0 + A(1,1) - A(2,2) - A(3,3)) * 2;
        q4 = (A(3,2) - A(2,3)) / S;
        q1 = 0.25 * S;
        q2 = (A(1,2) + A(2,1)) / S;  
        q3 = (A(1,3) + A(3,1)) / S; 
    elseif A(2,2) > A(3,3)
        S  = sqrt(1.0 + A(2,2) - A(1,1) - A(3,3)) * 2;
        q4 = (A(1,3) - A(3,1)) / S;
        q1 = (A(1,2) + A(2,1)) / S;
        q2 = 0.25 * S;  
        q3 = (A(2,3) + A(3,2)) / S; 
    else
        S  = sqrt(1.0 + A(3,3) - A(1,1) - A(2,2)) * 2;
        q4 = (A(2,1) - A(1,2)) / S;
        q1 = (A(1,3) + A(3,1)) / S;
        q2 = (A(2,3) + A(3,2)) / S;  
        q3 = 0.25 * S; 
    end


Q=[-q1 -q2 -q3 q4]';
Q_inv=[q1 q2 q3 q4]';
norm_q=q1^2+q2^2+q3^2+q4^2 %must be 1
x_axis=[1,0,0,0]';
y_axis=[0,1,0,0]';
z_axis=[0,0,1,0]';
CRF_x=qmul_i(qmul_i(Q,x_axis),Q_inv)
CRF_y=qmul_i(qmul_i(Q,y_axis),Q_inv)
CRF_z=qmul_i(qmul_i(Q,z_axis),Q_inv)

fprintf('Isat is %s \n\n',num2str((Isat*[1,1,1]')'));
fprintf('CoM is %s \n\n',num2str(CoM'));
fprintf('Q is %s \n\n',num2str(Q'));
fprintf('Mass of AAUSAT3: %f \n\n',mSat);

end
