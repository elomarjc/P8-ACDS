%%Albedo test
close all
clear
clc

fontsize=14;
n=20; %number of samples

figure
hold on

%plot n reflectivity matrixes in a contour plot
for i=0:(n-1)
    initial_time_JD=2455153.1954225+i*24/14.81823818/n/24;
    sim('albedo_test', 1);
    contour3(reflectivity(:,:,1),500)
    colormap jet
    %calculate total irradiance in percentage for the given position
    irradiance_total=0;
    for j=1:180
        for k=1:288
            irradiance_total=irradiance_total+reflectivity(j,k,1);
        end
    end
    irradiance(i+1)=irradiance_total/1366*100;
end
colorbar
irradiance

%simulate one orbit
sim('eclipse_test', 5830);
plot((sun_longitude+180)/1.25,sun_latitude+90,'.y');
h1=plot((moon_longitude+180)/1.25,moon_latitude+90,'.');
set(h1,'MarkerEdgeColor', [0.5 0.5 0.5]);
for i=1:length(eclipse)
    if(eclipse(i)==1)
        plot((sc_longitude(i)+180)/1.25,sc_latitude(i)+90,'.r')
    else
        plot((sc_longitude(i)+180)/1.25,sc_latitude(i)+90,'.b')
    end
    
end

%Add some labels
xlabel('X-axis cells [-]','FontSize',fontsize)
ylabel('Y-axis cells [-]','FontSize',fontsize)
zlabel('Cell irradiance [W/m^2]','FontSize',fontsize)

%Add a picture from orbitron with same orbit
rgb = imread('eclipse_map.jpg'); 
surface('XData',[0 288; 0 288],'YData',[0 0; 180 180],...
        'ZData',[0 0; 0 0],'CData',flipdim(rgb,1),...
        'FaceColor','texturemap','EdgeColor','none');

%Fix axis and view    
view(7,65);
grid on;
set(gcf, 'color', 'white');
xlim([0 288]);
ylim([0 180]);