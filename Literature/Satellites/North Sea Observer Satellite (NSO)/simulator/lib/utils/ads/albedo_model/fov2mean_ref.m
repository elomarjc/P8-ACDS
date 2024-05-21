% Calculates a new reflectivity matrix on the basis of 
% the albedo reflectivity matrix for the Earths surface.
% The new reflectivity matrix means all the reflectivity 
% points that are visible from the satellite on the basis 
% of the orbiting height.
%
% Param:
%   fovm: field of view matrix (from sat_fov_meridian)
%   or: The radius of the orbit from Earths center i meters. 
%       Must be larger that mean radius of earth (6371010m)
%   a_toms: is a TOMS measurement matrix 
%       180 x 288 reflectivity matrix loaded from epr_data
%   use_flux_mask: boolean flag - use flux or plain FOV as mask.
%
% Output:
%   mrm: reflectivity matrix (180 x 288) with data meaned 
%        on the basis of FOV
%   ones_count: number of measure points in FOV for each 
%        latitude along a meridian
%
function [mrm sum_matr ones_count] = fov2mean_ref(fovm,fov_flux,a_toms,or,use_flux_mask)
  global display_figs
  if ~exist('display_figs','var'),
    display_figs = true;
  end
  %conversion
  deg2rad = pi/180;
  rad2deg = 180/pi;
  
  if display_figs,
    %display for visual verification
    disp(['Displaying albedo FOV for visual verification (' datestr(now) ')...'])
    fig_h = figure(1);
    set(fig_h,'Name','Visual verification of mean computation')
    %subplot(2,1,1);
    img = reshape(a_toms,180,288);
    image(img*56);
    title('Albedo map')
    xlabel('longitude')
    ylabel('latitude')
    fov_fig_prop;%set common axes properties and stuff
    shg
    pause(.1)
  end
  
  status_disp = 12;% number of times to display a status line

  % allocating ram
  mrm = zeros(180,288);
  ones_count = zeros(180);

  %count nonzero entries
  disp(['Counting non-zero entries in FOV (' datestr(now) ')...'])
  fov_vec = reshape(fovm,180,180*288);
  ones_count = dot(fov_vec',fov_vec');
  
  %rotating toms 1--> 288 index
  disp([ 'Rotating TOMS data (' datestr(now) ')...'])
  a_rot_toms = zeros(288,180,288);
  for i_long_rot=1:288,
    %a_rot_toms(288+1-i_long_rot,:,:) = [a_toms(:,i_long_rot:288) a_toms(:,1:i_long_rot - 1)]; 
    a_rot_toms(i_long_rot,:,:) = [a_toms(:,i_long_rot:288) a_toms(:,1:i_long_rot - 1)]; 
  end
  
  %calculating sum (applying mask)
  disp([ 'Calculating sum in FOV (' datestr(now) ')...'])
  a_rot_toms_vec = reshape(a_rot_toms,288,180*288);
  fov_flux_vec =  reshape(fov_flux,180,180*288);

  %should flux be masked?
  if use_flux_mask,
    %using cos(rho) as mask
    sum_matr = fov_flux_vec*a_rot_toms_vec';
  else
    %using FOV as mask
    sum_matr = fov_vec*a_rot_toms_vec';
  end

  %computation of mean values
  disp(['Computing mean values (' datestr(now) ')...'])
  for i_long=1:288,
    mrm(:,i_long) = sum_matr(:,i_long)./ones_count';
  end
 
  
  
  
  
