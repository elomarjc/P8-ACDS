%Converts vector to sun in satellite bodyframe to angles of incident light
%for each sun sensor.
%
%  param: 
%   R_sun_sc : (1 x 3) vector : [x, y, z]
%
%  output:
%   (1 x 6) vector : [x+, y+, z+, x-, y-, z-]

function angles = sun_emu_vec2angles(R_sun_sc, name)
%Defining array with vectors normal to satellite planes
side = [1 0 0
        0 1 0 
        0 0 1];

%Output angles vector = x+,y+,z+,x-,y-,z-
angles = [90; 90; 90; 90; 90; 90];

%Calculating angle vector
sun_length = length(R_sun_sc);
sun_norm = norm(R_sun_sc);

%no input handling
if sun_norm == 0
    disp(['Sensor emulation error: Input vector to vec2angles from ' name ' has zero length! Returning']);
    return
end

for j = 1:sun_length,
  angle = 180/pi*(acos(norm(dot(R_sun_sc,side(:,j)))/(sun_norm)));
  if R_sun_sc(j) < 0   
    angles(j+3) = angle;
  else
    angles(j) = angle;
  end
end
 