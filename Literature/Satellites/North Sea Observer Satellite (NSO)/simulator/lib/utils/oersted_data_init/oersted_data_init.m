% This script is used for processing a 65 mb raw datafile from the Oersted
% satellite. The datafile contains samples approx every 1 sec. This is
% reduced by the factor step_size.
%
% For further description of format see:
% http://www.dmi.dk/projects/oersted/SDC/Format-MAG-L-1.0.shtml
% For further description of data see:
% http://www.dmi.dk/projects/oersted/SDC/DataDesc-2.3.shtml
% 
% From the datafile the following vectors is created:
%   Pos. data ECIRF
%   Pos. data ECEFRF
%   vel. data ECIRF
%   vel. data ECEFRF
%   mag. data ECIRF
%   mag. data ECEFRF
% Furthermore sample_number, the day and second of the measurement is also
% extracted
%
% All data is saved in the struct "oersted_data", which is saved to a file
% named "oersted_data.mat"
%
% This script uses a Simulink for data conversion. (It has dependencies!)
%% Preprocessing of Oersted datafile

% Format of datafile:
% Day,seconds,Radius,theta,phi,XXX,XXX,XXX,XXX,Mag_1_,Mag_2,Mag_3
% Where Radius,theta and phi is parameters for spherical coordinates

datafile=load('oersted-data/ml020210to17.txt');
stepsize = 5; 
% Init var. before the for-loop
OEdata_temp=zeros(round(length(datafile)/stepsize),6);
oersted_data.time.day=zeros(round(length(datafile)/stepsize),1);
oersted_data.time.second=zeros(round(length(datafile)/stepsize),1);
oersted_data.mag.ECEFRF=zeros(round(length(datafile)/stepsize),3);
datafile_marker_1=zeros(round(length(datafile)/stepsize),1);
datafile_marker_2=zeros(round(length(datafile)/stepsize),1);
sampletime=zeros(round(length(datafile)/stepsize),1);

stepsize=5;
day_tmp=datafile(1,1);
tmp3=1;
%find index for days
cnt=1;
for day_number=min(datafile(:,1)):1:max(datafile(:,1))
    day_index(cnt)=find(datafile(:,1)-day_number == 0,1);
    cnt=cnt+1;
end

%calculate mean Ørsted data sample time
cnt=1;
sample_rate_temp=zeros(40000,1);
for cnt=1:1:40000
    sample_rate_temp(cnt)=datafile(cnt+1,2)-datafile(cnt,2);
end
sample_rate_mean=mean(sample_rate_temp);

ii=0;
%Create time vector
for i=1200:stepsize:length(datafile)-stepsize,
    ii=ii+1;
    	oersted_data.time.day(ii)=day_tmp+fix(i/86400);
    oersted_data.time.second(ii)=mod(i,(oersted_data.time.day(ii)-day_tmp)*86400);
end

ii=1;
for i=1200:stepsize:length(datafile)-stepsize,
    ii=ii+1;

    
     tmp_1=abs(round(i/sample_rate_mean)-1000);
     datafile_marker_1(ii) = tmp_1 - 1 + find(datafile(tmp_1:tmp_1+2000,2)-sec(ii) > 0 , 1);
     datafile_marker_2(ii) = tmp_1 - 1 + find(datafile(tmp_1:tmp_1+2000,2)-sec(ii+1) > 0 , 1);
    

 
    if datafile_marker_1(ii)==datafile_marker_2(ii) 
        datafile_marker_2(ii)=datafile_marker_1(ii)+1;
    end

     if datafile_marker_1(ii)==datafile_marker_1(ii-1) %check timestep in datafile and correct stepsize
         tmp3=tmp3+1;
     else
         tmp3=1;
     end

    sampletime(ii) = (datafile(datafile_marker_2(ii),2)-datafile(datafile_marker_1(ii),2));
    if sampletime(ii) < 0     % Check for day shift
        ii
        sampletime(ii) = (datafile(datafile_marker_2(ii),2) + 86400 - datafile(datafile_marker_1(ii),2));
    end

    OEdata_temp(ii,1)=datafile(datafile_marker_1(ii),3)+((datafile(datafile_marker_2(ii),3)-datafile(datafile_marker_1(ii),3))/sampletime(ii))*stepsize*tmp3;
    OEdata_temp(ii,2)=datafile(datafile_marker_1(ii),4)+((datafile(datafile_marker_2(ii),4)-datafile(datafile_marker_1(ii),4))/sampletime(ii))*stepsize*tmp3;
    OEdata_temp(ii,3)=datafile(datafile_marker_1(ii),5)+((datafile(datafile_marker_2(ii),5)-datafile(datafile_marker_1(ii),5))/sampletime(ii))*stepsize*tmp3;
    OEdata_temp(ii,4)=datafile(datafile_marker_1(ii),10)+((datafile(datafile_marker_2(ii),10)-datafile(datafile_marker_1(ii),10))/sampletime(ii))*stepsize*tmp3;
    OEdata_temp(ii,5)=datafile(datafile_marker_1(ii),11)+((datafile(datafile_marker_2(ii),11)-datafile(datafile_marker_1(ii),11))/sampletime(ii))*stepsize*tmp3;
    OEdata_temp(ii,6)=datafile(datafile_marker_1(ii),12)+((datafile(datafile_marker_2(ii),12)-datafile(datafile_marker_1(ii),12))/sampletime(ii))*stepsize*tmp3;
    
end


i=1;
%% old version
% for ii=1:stepsize:length(datafile),
%     %sperical pos. data:
%     OEdata_temp_old(i,1)=datafile(ii,3);
% 	OEdata_temp_old(i,2)=datafile(ii,4);
% 	OEdata_temp_old(i,3)=datafile(ii,5);
%     
%     %sperical mag. data
%     OEdata_temp_old(i,4)=datafile(ii,10);
% 	OEdata_temp_old(i,5)=datafile(ii,11);
% 	OEdata_temp_old(i,6)=datafile(ii,12);
%     
% 	oersted_data.time.day(i)=datafile(ii,1);
%     oersted_data.time.second(i)=datafile(ii,2);
%    
%     i=i+1;
% end
Put sample info in struct. (needed by Simulink)
    oersted_data.time.sample_number=[1:length(OEdata_temp)]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Convert Sperical pos. data to cartesian coordinates and shift reference frame from ECEFRF to ECIRF for mag. data and pos. data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For the ease, a Simulink file is used. The result is put in OEdata_temp3
%sim('oersted_data_init_sim', [1 length(OEdata_temp)]);
sim('oersted_data_init_sim', [1 50000]);

% put the results in the struct
oersted_data.pos.ECEFRF=zeros(length(OEdata_temp3),3);
oersted_data.pos.ECIRF=zeros(length(OEdata_temp3),3);
oersted_data.mag.ECIRF=zeros(length(OEdata_temp3),3);
oersted_data.mag.ECEFRF=zeros(length(OEdata_temp3),3);

for i=1:1:length(OEdata_temp3),
    oersted_data.pos.ECEFRF(i,1)=OEdata_temp2(i,1);
    oersted_data.pos.ECEFRF(i,2)=OEdata_temp2(i,2);
    oersted_data.pos.ECEFRF(i,3)=OEdata_temp2(i,3);

    oersted_data.pos.ECIRF(i,1)=OEdata_temp3(i,1);
    oersted_data.pos.ECIRF(i,2)=OEdata_temp3(i,2);
    oersted_data.pos.ECIRF(i,3)=OEdata_temp3(i,3);

    oersted_data.mag.ECEFRF(i,1)=OEdata_temp4(i,1);    
    oersted_data.mag.ECEFRF(i,2)=OEdata_temp4(i,2);    
    oersted_data.mag.ECEFRF(i,3)=OEdata_temp4(i,3);    

    oersted_data.mag.ECIRF(i,1)=OEdata_temp5(i,1);    
    oersted_data.mag.ECIRF(i,2)=OEdata_temp5(i,2);    
    oersted_data.mag.ECIRF(i,3)=OEdata_temp5(i,3);    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculation of velocity %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


oersted_data.vel.ECIRF = zeros(length(oersted_data.pos.ECIRF),3);
sample_time_temp = zeros(length(oersted_data.pos.ECIRF),1);

for i=1:1:length(oersted_data.pos.ECIRF),

    sampletime = (oersted_data.time.second(i+1)-oersted_data.time.second(i));
    if sampletime < 0     % Check for day shift
        sampletime = (oersted_data.time.second(i+1) + 86400 - oersted_data.time.second(i));
    end
    sample_time_temp(i)=sampletime;
    oersted_data.vel.ECIRF(i,1) = (oersted_data.pos.ECIRF(i+1,1) - oersted_data.pos.ECIRF(i,1)) / sampletime;
    oersted_data.vel.ECIRF(i,2) = (oersted_data.pos.ECIRF(i+1,2) - oersted_data.pos.ECIRF(i,2)) / sampletime;
    oersted_data.vel.ECIRF(i,3) = (oersted_data.pos.ECIRF(i+1,3) - oersted_data.pos.ECIRF(i,3)) / sampletime;

    if i+1 == length(oersted_data.pos.ECIRF)
        oersted_data.vel.ECIRF(i+1,1) = oersted_data.vel.ECIRF(i,1);
        oersted_data.vel.ECIRF(i+1,2) = oersted_data.vel.ECIRF(i,2);
        oersted_data.vel.ECIRF(i+1,3) = oersted_data.vel.ECIRF(i,3);
        break;
    end

end

%% Clear variables and save the struct
clear datafile;
clear stepsize;
clear i;
clear ii;
clear sample_time_temp;
clear sampletime;
clear tout;
clear OEdata_temp;
clear OEdata_temp2;
clear OEdata_temp3;
clear OEdata_temp4;
clear OEdata_temp5;

save ../oersted_data oersted_data