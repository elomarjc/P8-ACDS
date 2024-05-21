%Used for dimensioning of the camera system to set requirements

resolution = [1240 1024];
oil_size =150; %quadratic area in meters 
oil_res = [4 4]; %required size of oil spill is pixel 
orbit_height = 700e3; %m

max_pixel_move_oil = min(resolution)/2-max(oil_res)/2;

ground_size = max_pixel_move_oil/max(oil_res)*oil_size;

max_dev_ang = rad2deg(atan(ground_size/orbit_height))