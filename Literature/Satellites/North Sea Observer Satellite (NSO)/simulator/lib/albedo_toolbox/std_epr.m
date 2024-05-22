% STD_EPR Calculate standard deviation of epr data in MAT files..
%
% [std_val, std_lat] = std_epr(files, mean_val, mean_lat)
%
% files is a string of a wildcard filename or a directory. Only files with
% '.mat' extension are read. mean_val and mean_lat are the mean values of
% the epr data in files. std_val contains the standard deviation values
% and std_lat contains the standard deviation values over latitudes.
% std_val.start_time is the earliest file start time, and std_val.stop_time
% is the latest. std_val is an epr struct. std_lat is a vector.
%
% $Id: std_epr.m,v 1.1 2006/03/17 12:13:52 mnkr02 Exp $

function [std_val, std_lat] = std_epr(file, mean_val, mean_lat)

%%%%%%%%%%%%%%%%%%
% Init
%%%%%%%%%%%%%%%%%%

[parm_path,name,ext,versn] = fileparts(file);
files = dir([file]);
nfiles = size(files,1);

if isdir(file)
  parm_path = file;
end

% Allocate accumulated variance array
epr_acc = zeros(180,288);

% Allocate accumulated variance array with 180x1 mean array for simplified
% latitude based epr model
epr_lat_acc = zeros(180,288);

% Allocate valid sample count
epr_count = zeros(180,288);

% Initialize start and stop times
start_time = zeros(nfiles,1);
stop_time = zeros(nfiles,1);

%%%%%%%%%%%%%%%%%%
% Calculate standard deviation
%%%%%%%%%%%%%%%%%%

fprintf(1,'Progress: %3.0f%%',0);

for ifile = 1:nfiles
  [pathstr,name,ext,vers] = fileparts(files(ifile).name);
  if strcmp(ext,'.mat')
    % Read next file
    epr = load(fullfile(parm_path, files(ifile).name));
    % Get start and stop times
    start_time(ifile) = epr.start_time;
    stop_time(ifile) = epr.stop_time;
    for index = 1:180*288
      if ~isnan(epr.data(index))
        % Accumelate data values if valid
        epr_acc(index) = epr_acc(index) + (epr.data(index)-mean_val.data(index))^2;
        epr_lat_acc(index) = epr_lat_acc(index) + (epr.data(index)-mean_lat(mod(index-1,180)+1))^2;
        % Track valid sample count for mean value calculation
        epr_count(index) = epr_count(index) + 1;
      end
    end
    start_time(ifile) = inf;
  end
  fprintf(1,'\b\b\b\b%3.0f%%',ifile*100/nfiles);
end

% Return std_val
warning('off');
std_val = epr_struct(sqrt(epr_acc./epr_count),min(start_time),max(stop_time),'Standard Deviation');
if min(min(epr_count)) == 0
  fprintf(1,'\nstd_val contains %.0f%% empty values (NaN).',sum(sum(isnan(std_val.data)))*100/(180*288));
end

% Return std_lat
if nargout > 1
  std_lat = sqrt(sum(epr_lat_acc')./sum(epr_count'))';
  if min(sum(epr_count')) == 0
    fprintf(1,'\nstd_lat contains %.0f%% empty values (NaN).',sum(isnan(std_lat))*100/180);
  end
end

fprintf('\n');