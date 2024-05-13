%% Eclipse, sun and moon test
close all
clear
clc

fontsize=12;

%simulate one orbit
sim('eclipse_test', 5830);

%plot map showing start and satellite trajectory
h = subplot(2,1,1);
set(h,'position',[0.1 0.53 0.85 0.45]); 
hold on
rgb = imread('eclipse_map.jpg'); %picture generated with orbitron
image([-180 180],[90 -90],rgb);
set(gca,'YDir','normal')
set(gca,'ytick',[-90 -60 -30 0 30 60 90])
set(gca,'xtick',[])
xlim([-180 180]);
ylim([-90 90]);
for i=1:length(eclipse)
    if(eclipse(i)==1)
        plot(sc_longitude(i),sc_latitude(i),'.r')
    else
        plot(sc_longitude(i),sc_latitude(i),'.b')
    end
    
end
plot(sun_longitude,sun_latitude,'.y');
h1=plot(moon_longitude,moon_latitude,'.');
set(h1,'MarkerEdgeColor', [0.5 0.5 0.5]);
ylabel('Latitude [deg]','FontSize',fontsize);

%plot map showing end and satellite trajectory
h=subplot(2,1,2);
set(h,'position',[0.1 0.05 0.85 0.45]); 
hold on
rgb = imread('eclipse_map2.jpg'); %picture generated with orbitron
image([-180 180],[90 -90],rgb);
set(gca,'YDir','normal')
set(gca,'ytick',[-90 -60 -30 0 30 60 90])
set(gca,'xtick',[-180 -150 -120 -90 -60 -30 0 30 60 90 120 150 180])
xlim([-180 180]);
ylim([-90 90]);
for i=1:length(eclipse)
    if(eclipse(i)==1)
        plot(sc_longitude(i),sc_latitude(i),'.r')
    else
        plot(sc_longitude(i),sc_latitude(i),'.b')
    end
   
end
plot(sun_longitude,sun_latitude,'.y');
h1=plot(moon_longitude,moon_latitude,'.');
set(h1,'MarkerEdgeColor', [0.5 0.5 0.5]);
xlabel('Longitude [deg]','FontSize',fontsize);
ylabel('Latitude [deg]','FontSize',fontsize);

i=0;
for j=1:length(eclipse)
    if(eclipse(j,1)==1)
        i=i+1;
    end
end
i/60