%%Albedo total irradiance plot for one orbit
close all
clear
clc

fontsize=14;
n=200; %number of samples

figure
hold on

%plot irradiance in percent
for i=0:(n-1)
    initial_time_JD=2455153.1954225+i*24/14.81823818/n/24;
    sim('albedo_test', 1);
    %calculate total irradiance in percentage for the given position
    irradiance_total=0;
    for j=1:180
        for k=1:288
            irradiance_total=irradiance_total+reflectivity(j,k,1);
        end
    end
    irradiance=irradiance_total/1366*100;
    plot3([sc_longitude(1) sc_longitude(1)],[sc_latitude(1) sc_latitude(1)],[0 irradiance],'-','LineWidth',2,'Color',[.5 .5 .5])
end

%simulate one orbit
sim('eclipse_test', 5830);
plot(sun_longitude,sun_latitude,'.y');
h1=plot(moon_longitude,moon_latitude,'.');
set(h1,'MarkerEdgeColor', [0.5 0.5 0.5]);
for i=1:length(eclipse)
    if(eclipse(i)==1)
        plot(sc_longitude(i),sc_latitude(i),'.r')
    else
        plot(sc_longitude(i),sc_latitude(i),'.b')
    end
    
end

%Add some labels
xlabel('Longitude [deg]','FontSize',fontsize)
ylabel('Latitude [deg]','FontSize',fontsize)
zlabel('Reflected irradiance [%]','FontSize',fontsize)

%Add a picture from orbitron with same orbit
rgb = imread('eclipse_map.jpg'); 
surface('XData',[-180 180; -180 180],'YData',[-90 -90; 90 90],...
        'ZData',[0 0; 0 0],'CData',flipdim(rgb,1),...
        'FaceColor','texturemap','EdgeColor','none');

%Fix axis and view    
view(7,65);
grid on;
set(gcf, 'color', 'white');
xlim([-180 180]);
ylim([-90 90]);
zlim([0 50]);