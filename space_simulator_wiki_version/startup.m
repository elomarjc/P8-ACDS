function startup(mex)
%this script is automatically run when matlab is started in this dir
% it adds the following to the matlab path
sub_dirs = {
    'lib/'
    'lib/AlbedoToolbox-1.0'
    'lib/AlbedoToolbox-1.0/refl_data'
    'lib/AlbedoToolbox-1.0/refl_data/2005'
    'lib/utils'
    'lib/utils/ephemeris'
    'lib/utils/math_utils'
    'lib/utils/animation'
    'lib/utils/ADS'
    'lib/utils/ADS/AAUSAT-II'
    'lib/utils/ADS/AAUSAT3'
    'lib/utils/ADS/SATBALL'
    'lib/pictures'
    'misc/'
    'design/'
    'models/'
    'test/'
    'test/IGRF/'
    'test/SGP4/'
    'test/albedo_eclipse/'
    'test/ads/'
    'test/b_dot/'
    'test/Bdot_control/'
    'test/animation/'
%     'test/animation/movie/'
    'test/disturbance/'
    'test/sm_control'
    'test/satball'
	   };

len=length(sub_dirs);
disp('Adding to path...')
for i= 1:len,
  temp = [pwd '/' char(sub_dirs(i))];
% disp(temp)
  addpath(temp)
end
disp('Done!')

clear sub_dirs len i temp
if nargin~=0
disp('Compiling mex files...')
cd misc;
compile_mex;
cd ..;
disp('Done!')
end

