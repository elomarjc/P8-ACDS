function startup(mex)
%this script is automatically run when matlab is started in this dir
% it adds the following to the matlab path
sub_dirs = {
    'lib/'
    'lib/utils'
    'lib/utils/ads'
    'lib/utils/ads/albedo_model'
    'lib/utils/acs'
    'lib/utils/ephemeris'
    'lib/utils/sensor_emulation'
    'lib/utils/math_utils'
    'lib/utils/animation'
    'lib/pictures'
    'lib/albedo_toolbox'
    'lib/albedo_toolbox/epr_data/'
    'lib/albedo_toolbox/epr_data/2001'
    'lib/albedo_toolbox/epr_data/2002'
    'lib/albedo_toolbox/epr_data/2003'
    'misc/'
    'design/'
    'models'
    'test'
    'test/m-files'
    ''
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

