% Some stupid FOV calculations for on-board usage

% Earth radius
rE = 6371.01e3;

% Altitude
alt = 650e3;

% Discretize resolution (degrees)
dis = 2;

% FOV
fov = round(90 - asin(rE/(rE+alt))*180/pi);

% Generate FOV area table
steps = abs(2*floor(fov / dis));
f = -1*fov;
fovtab(1)=0;
xax(1)=f;
for i = 0:steps
    f=f+dis;
    fovtab(i+2) = cos(f/fov * pi/2);
    xax(i+2)=f;
end

% Plot it
figure(1)
stem(xax,fovtab);
grid;
xlabel('Current FOV lattitude [deg]');
ylabel('Area compensation [eu]');

% Let's combine the albedo and FOV table
load albedotab;

% It's quite easy - convolution and mean the result
obalbedotab = convn(ads_c_mean, fovtab, 'same');
obalbedotab = obalbedotab/steps;

% Plot it
figure(2)
stem((-45:44)*2,obalbedotab*100);
grid;
xlabel('Lattitude [deg]');
ylabel('FOV corrected reflectivity % [eu]');
