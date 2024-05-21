% Calculates FOV matrix from satellite position
%
% Parameters:
%   long_0: longitude of the satellite position (radians)
%   lat_0:  lattitude of the satellite position (radians)
%   or: The radius of the orbit from Earths center i meters. 
%       Must be larger that mean radius of earth (6371010m) 
%
% Output: An FOV matrix with the dimentions of (180 x 288)
%         Entries are '1' for visible points and
%         '0' for non-visible points
%
% to test run:
%
%  long_0 = 0, lat_0=0, or=6371010+5e5, long=deg2rad([-179.375:1.25:179.375]); 
%  sat_fov(long_0, lat_0, long, or);
function [fovm fov_flux] = sat_fov(long_0, lat_0, long, or)
  persistent r_earth arc_theta lat_fov_index;
  
  %allocation of memory
  fovm = zeros(180,288);
  fov_flux = fovm;
  
  %conversion
  deg2rad = pi/180;
  rad2deg = 180/pi;
    
  %Radius of earth in meters.
  r_earth = 6371010; 
  
  %Angle of visible arc measured from earths center 
  %alpha_limit = acos(r_earth/or);
  cos_rho_0 = r_earth/or;%cos(alpha_limit);

  r_o = or;
  r_e = r_earth;
  
  % new!!!!
  lat_fov_index=1;
  for lat = deg2rad*(-89.5:1:89.5),
    % small circle equation (Spheric section)
    cos_rho = sin(lat)*sin(lat_0) + cos(lat)*cos(lat_0)*cos(long-long_0);
    
    % fields of view mask for input satellite vector
    % mask is 1 inside circle, 0 elsewhere
    fovm(lat_fov_index, :) = double(cos_rho_0 < cos_rho);
    
    %field of view flux mask
    cos_beta = (r_o - r_e*cos_rho) /(sqrt(r_e^2+r_o^2 - 2*r_e*r_o*cos_rho));
    
    cos_theta = cos(acos(cos_beta) + acos(cos_rho));
    
    if max(cos_rho) > 1,
      error(['cos(rho) is greater than zero, cannot proceed. Latitude is ' num2str(rad2deg*(lat)) ]);
    end
    
    %save flux masked with fov
    fov_flux(lat_fov_index, :) = cos(acos(cos_beta)+ acos(cos_rho)) .* fovm(lat_fov_index, :);
    
    %increment interation
    lat_fov_index = lat_fov_index + 1;
  end
  
  

% $$$   % old!!!!
% $$$   lat_fov_index=1;
% $$$   for lat = deg2rad*(-89.5:1:89.5),
% $$$     %fields of view mask for input satellite vector
% $$$     fovm(lat_fov_index, :) = double(cos_rho_0 < (sin(lat)*sin(lat_0) + cos(lat)*cos(lat_0)*cos(long-long_0)));
% $$$     lat_fov_index = lat_fov_index + 1;
% $$$   end
% $$$  