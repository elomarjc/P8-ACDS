%Calculates 4'th order polinomial fit from measured data in file: 'sundata.csv'
%
%  param: none
%
%  output:
%   (1 x 5) polynomial coefficient vector for use with polyval.
%
% see also: polyfit, polyval

function p = sun_emu_model()
%reads data into aray from file
% column names:
%VPOSACTH_EU,HPOSACTH_EU,CUR1_EU,CUR2_EU,CUR3_EU,TEMP1_EU,TEMP2_EU,TEMP3_EU
[sundata] = textread('sun_emu_data.csv','','delimiter',',','headerlines',1);

%copy data to new variables
angle = sundata(:,2);   %inclination angles
measurement = sundata(:,3); %data from one sun sensor

j = 1;
k = 1;
for f = 0:10:90,
    temp = 0;
    l = 0;
    while angle(k) == f,
        temp = temp + measurement(k)/1000;
        if k < 541,
            k = k+1;
        else
            break;
        end
        l = l+1;
    end
    sun_mean(j,1) = f;
    sun_mean(j,2) = temp/l;
    sun_norm(j,1) = f;
    sun_norm(j,2)=sun_mean(j,2)/sun_mean(1,2);
    j = j+1;
end
[p w] = polyfit(sun_mean(:,1),sun_mean(:,2),4);

