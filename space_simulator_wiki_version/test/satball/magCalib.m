% Clear command window
clc

% Close all figures
close all

% Save magnetometer data to 'magDataRaw.mat'
save('magDataRaw.mat','sat_mag');

% Extract magnetometer data along x, y, and z axes
mag_x = sat_mag.Data(:,1);
mag_y = sat_mag.Data(:,2);
mag_z = sat_mag.Data(:,3);

% Plot magnetometer data in 3D
figure;
plot3(mag_x,mag_y,mag_z,'*');
axis equal

% Combine magnetometer data into a matrix H
H = [mag_x mag_y mag_z];
% Display the size of matrix H
size(H)

% Perform magnetometer calibration
[U, c] = MgnCalibration(H);

% Correct magnetometer readings using calibration parameters
for i = 1:size(H,1)
    w(i,:) = (U*(H(i,:)' - c))';
end

% Plot calibrated magnetometer data in 3D
figure;
plot3(w(:,1),w(:,2),w(:,3),'*');
axis equal

% Set the vector v_mag
v_mag = [0.1145 0.3459 0.8283];
