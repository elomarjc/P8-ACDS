% Displays field of view from satellite 
% at 180 latitudes sequentially
%
% Params:
%    fovm: field of view matrix (i x 180 x 288)
%          wher i is arbitrary dimension
%
%
function fov_disp(fovm)

  %input checking
  img_count = size(fovm,1);
  if (ndims(fovm) ~= 3) || logical(sum(size(fovm) ~= [img_count 180 288])),
    disp(['fovm has wrong dimension, must be of size (i x 180 x 288). '...
	  'aborting display of image'])
    return
  end

  %displays the image
  figure(1);
  shg
  for lat_fov_index= 1:5:img_count,
    image(uint8(reshape(fovm(lat_fov_index,:,:),180,288)*56));
    title('Field of view from satellite')
    xlabel('longitude')
    ylabel('latitude')
    %set(gca,'DataAspectRatio',[1 1 1])
    set(gca,'xtick',[1  72 144 216 288],'XTickLabel',{'-180','-90','0','90','180'})
    set(gca,'ytick',[1 45 90 135 180],'YTickLabel',{'-90','-45','0','45','90'})
    pause(.01);
  end
