% Calculates FOV matrices along the meridian described by long_0
%
% Parameters:
%   long_0: longitude of the satellite position (radians)
%   or: The radius of the orbit from Earths center i meters. 
%       Must be larger that mean radius of earth (6371010m) 
%
% Output: An FOV matrix with the dimentions of (180 x 180 x 288)
%         Entries are '1' for visible points and
%         '0' for non-visible points
%
function [fovm fov_flux] = sat_fov_meridian(long_0,or)

  if exist('display_figs','var'),
     global display_figs;
  else
    disp('Diplay of figures is on (sat_fov_meridian.m)')
    display_figs = true;
  end
  
  %conversion
  deg2rad = pi/180;
  rad2deg = 180/pi;

  
  %struct size
  m_size = [180,180,288];
  
  %calculate field of view along meridian: long_0, (lat: -90 --> 90)
  disp(['Calculating field of view (' datestr(now) ')...'])
  long = deg2rad*(-179.375:1.25:179.375);
  long_0=-179.375;
  fovm = zeros(m_size);
  lat_fov_index = 1;
  for lat_0 = deg2rad*(-89.5:1:89.5),
    if ~mod(lat_fov_index,20)
      disp([ '... processing latitude index: ' ... 
	     num2str(lat_fov_index) ' (' datestr(now) ')'])
    end
    [fovm_lat fov_flux_lat] = sat_fov(long_0, lat_0, long, or);
    fovm(lat_fov_index,:,:) = fovm_lat;
    fov_flux(lat_fov_index,:,:) = fov_flux_lat;
    lat_fov_index = lat_fov_index + 1;
  end
  disp(['...fov calculation done (' datestr(now) ')'])

  
  if display_figs,
    %displaying field of view
    disp(['Displaying field of view ' datestr(now) ')'])
    figure;
    subplot(2,1,1);
    fov_disp(fovm);
    title('FOV mask')
    subplot(2,1,1);
    fov_disp(fov_flux);
    title('FOV flux distribution')
    disp(['done (' datestr(now) ')'])
  end
