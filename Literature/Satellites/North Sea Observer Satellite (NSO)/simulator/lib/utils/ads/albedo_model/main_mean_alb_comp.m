%script to process data and generate .cvs files containing meaned data

%toggle debug modes
if ~exist('force_recalc','var'),
  disp('setting up debug modes')
  display_figs = false
  force_recalc = false
  force_reload = false
  use_bogus = false
end


time_start = now;
disp(['Simulation processing(' datestr(time_start) ')...'])
%disp(['Debugging mode? ' logical2string(debug)])
%disp(['Forcing recalculation of all variables? ' logical2string(force_recalc)])

%load data
if exist('fov_sim_all.mat','file'),
  if ~force_recalc ,
    disp(['Loading data from previous calculation (' datestr(time_start) ')...'])
    load fov_sim_all
  end
end
output_updated = false;

%conversion
deg2rad = pi/180;
rad2deg = 180/pi;

%FOV meridian - MUST MATCH TOMS data!!!
fov_meridian_0 = -179.375;

if ~exist('or','var'),
  %orbit radius
  or = 6371010+5e5;
end

% create FOV matrix
if force_recalc || ~exist('fovm','var') || ~exist('fov_flux','var'),
  if ~(force_recalc ...
       || ~exist('fov_sim_fovm.mat','file') ...
       || ~exist('fov_sim_fovm.mat','file')),
    disp(['Loading FOV data from file (' datestr(now) ')...'])
    load fov_sim_fovm
    load fov_sim_fov_flux
  else
    disp(['Recalculating FOV along meridian (' datestr(now) ')...'])
    [fovm fov_flux]=sat_fov_meridian(fov_meridian_0,or);
    save fov_sim_fovm fovm
    save fov_sim_fov_flux fov_flux
    output_updated = true;
  end
  if ~(exist('fovm','var') == 1),
    error('Could not load FOV data, exiting')
  end
end




if display_figs,
  %displaying field of view
  disp(['Displaying field of view (' datestr(now) ')...'])
  fig_h=figure(1);
  set(fig_h,'Name','Field of view from satellite')
  pause(.01);
  for lat_fov_index= 1:5:180,
    subplot(2,1,1);
    %fov_disp(fovm);
    image(reshape(fovm(lat_fov_index,:,:),180,288)*56);
    title('Field of view from satellite')
    fov_fig_prop;
    subplot(2,1,2);
    image(reshape(fov_flux(lat_fov_index,:,:),180,288)*56);
    %fov_disp(fov_flux);
    title('FOV flux distribution')
    fov_fig_prop;
    pause(.01);
  end%%%%%%%%%%%%%%
end



 
	     

if force_recalc || force_reload || ~exist('a_toms','var') ,
  %loading Earth reflectivity data from TOMS experiment
  if ~use_bogus,
    disp(['Loading Earth reflectivity data (TOMS format) (' datestr(now) ')...'])
    a_toms = load('../../../space_environment_emulation/albedo/toolbox/epr_data/2003/ga030101-031231.mat');
    a_toms = a_toms.data;
  else
    %test data:
    disp(['Generating bogus reflectivity data (TOMS format) (' datestr(now) ')...'])
    a_toms = zeros(180,288);
    for i=1:90, 
      a_toms(i,i) = i;
    end
  end
  output_updated = true;
end

if force_recalc || ~(exist('mrm','var') && exist('ones_count','var')),
  %calculate mean data
  disp(['Calculating meaned data (' datestr(now) ')...'])
  if ~exist('use_flux_mask','var'),
    use_flux_mask = true;
  end
  [mrm sum_matr ones_count] = fov2mean_ref(fovm,fov_flux,a_toms,or, use_flux_mask);
  disp(['...mean calculation done (' datestr(now) ')'])
  output_updated = true;
  if exist('mrm','var') && exist('ones_count','var'),
    %save results on succes
    file = ['meaned_ref_matrix_' num2str(or) '.csv'];
    csvwrite(file, mrm);
    save fov_sim_mrm mrm
    disp(['Saved meaned albedo matrix (mrm) in : ' file ', and fov_sim_mrm.mat'])
    file = ['measure_mask_size_' num2str(or) '.csv'];
    csvwrite(file, ones_count);
    save fov_sim_ones_count ones_count
    disp(['Wrote FOV count matrix (ones_count) : ' file ', and fov_sim_ones_count.mat'])
  else
      disp('No result from simulation :-(')
  end
end

%display stats
disp(['Simulation started ' datestr(time_start) ', ended ' datestr(now)]);

if output_updated || ~exist('fov_sim_all.mat','file')
  if (exist('or','var') && exist('mrm','var') && exist('ones_count','var') && exist('a_toms','var')),
    save fov_sim_all or a_toms fovm  ones_count mrm 
    disp('Save succesful :-)')
  else
    error('No sufficient data for saving')
  end
end

if display_figs,
  %displaying mean reflectivity
  disp(['Displaying results (' datestr(now) ')...'])
  fig_h = figure(2);
  px_dim=2;
  py_dim=3;
  cam_pos = 1.0e+03*[-1.2197 -1.0900 0.0017]+[-2005 100 -.2]
  shg
  set(fig_h,'Name','FOV Results (albedo, fov, mrm picture, mrm graph)')

  %albedo
  subplot(py_dim,px_dim,1);
  image(a_toms*56);
  title('Albedo input map')
  xlabel('longitude')
  ylabel('latitude')
  fov_fig_prop;%set common axes properties and stuff

  %illustrate procedure   
  for i_long_rot=1:36:288,
    for i_lat=1:10:180,
      if ~mod(i_lat-1,30)
	disp([ '... rotating displays'])
	
	% albedo mesh
	subplot(py_dim,px_dim,2);
	meshc(a_toms);
	set(gca,'CameraPosition', cam_pos);
	title('Albedo input 3D map')
	xlabel('longitude')
	ylabel('latitude')
	fov_fig_prop;%set common axes properties and stuff

	%Meaned reflectivity picture
	subplot(py_dim,px_dim,6)
	meshc(mrm); 
	set(gca,'CameraPosition', cam_pos);
	title('Meaned reflectivity')
	fov_fig_prop;

	%Summed reflectivity
	subplot(py_dim,px_dim,4);
	meshc(sum_matr);
	xlabel('Latitude [\circ]')
	ylabel('Reflectivity [%]')
	set(gca,'CameraPosition', cam_pos);
	pause(.01);
      end
      cam_pos = cam_pos + [-100 0 0];
      % 1 rotate background (fov2mean_ref: line 59)
      imgbg = [a_toms(:,i_long_rot:288) a_toms(:,1:i_long_rot - 1)];
      % 2 aply mask  (fov2mean_ref: line 66)
      img_mask = imgbg.* (reshape(fovm(i_lat,:,:),180,288)+.25)*0.8;
      img_flux = imgbg.* (reshape(fov_flux(i_lat,:,:),180,288)+.25)*0.8;
      % 3 rotate back (for viewing consistency)
      i=i_long_rot;
      n=288;
      %cant be done backwards, so first we mirror...
      img_mask = [img_mask(:,n:-1:1)];
      img_flux = [img_flux(:,n:-1:1)];
      %roll forward...
      img_mask = [img_mask(:,i:n) img_mask(:,1:i-1)];
      img_flux = [img_flux(:,i:n) img_flux(:,1:i-1)];
      %and mirror again
      img_mask = [img_mask(:,n:-1:1)];
      img_flux = [img_flux(:,n:-1:1)];

      %display procedure using flux...
      subplot(py_dim,px_dim,3);
      image(img_mask*56);
      title('Albedo FOV using fov mask')
      fov_fig_prop;
      
      %display procedure using flux...
      subplot(py_dim,px_dim,5);
      image(img_flux*56);
      title('Albedo FOV using flux compensated mask')
      fov_fig_prop;
      pause(.5);
    end
  end
  disp(['done (' datestr(now) ')'])
end
 
