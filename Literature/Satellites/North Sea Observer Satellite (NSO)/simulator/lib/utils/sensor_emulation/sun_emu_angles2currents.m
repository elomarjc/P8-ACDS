%Converts angles of incident light for each sun sensor 
%to the current outputs of the corresponding sensor.
%(min. value of currents are the sensor dark current)
%
%  param: 
%   angles : (1 x 6) vector : [x+, y+, z+, x-, y-, z-]
%   angle2current_model : model polynomial coefficient vector
%   i_dark : scalar value of the dark current in the sensor model
%
%  output:
%   (1 x 6) vector : [I_x+, I_y+, I_z+, I_x-, I_y-, I_z-]
function currents = sun_emu_angles2currents(angles, angle2current_model, i_dark)
currents = [0;0;0;0;0;0];
angle_length = length(angles);
for j = 1:angle_length,
    temp = abs(angles(j));
    if temp >= 90
        currents(j) = i_dark;
    else
        %%%%%%%%%%%%%%%%%%
        %% Cosine model %%
        %%%%%%%%%%%%%%%%%%
        %currents(j) = cos(deg2rad(temp));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Polynomial model based on Ørsted sun sensor data %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        currents(j) = polyval(angle2current_model, temp); 
    end
end

 