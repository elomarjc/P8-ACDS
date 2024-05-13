% MEAN_EPR Calculate mean of epr data in MAT files..
%
% [mean_val, mean_lat] = mean_epr(files)
%
% files is a string of a wildcard filename or a directory. Only files with
% '.mat' extension are read. mean_val contains the mean values and mean_lat
% contains the mean values meaned over latitudes. mean_val.start_time is the
% earliest file start time, and mean_val.stop_time is the latest. mean_val
% is an epr struct. mean_lat is a vector.
%
% $Id: mean_epr.m,v 1.1 2006/03/17 12:13:51 mnkr02 Exp $

function [mean_val, mean_lat] = mean_epr(file)

%%%%%%%%%%%%%%%%%%
% Init
%%%%%%%%%%%%%%%%%%

[parm_path,name,ext,versn] = fileparts(file);
files = dir([file]);
nfiles = size(files,1);

if isdir(file)
  parm_path = file;
end

% Allocate valid sample count
epr_acc = zeros(180,288);

% Allocate mean value array
epr_count = zeros(180,288);

% Initialize start and stop times
start_time = zeros(nfiles,1);
stop_time = zeros(nfiles,1);

%%%%%%%%%%%%%%%%%%
% Calculate mean
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
    % Get epr data
    for index = 1:180*288
      if ~isnan(epr.data(index))
        % Accumelate data values if valid
        epr_acc(index) = epr_acc(index) + epr.data(index);
        % Track valid sample count for mean value calculation
        epr_count(index) = epr_count(index) + 1;
      end
    end
  else
    start_time(ifile) = inf;
  end
  fprintf(1,'\b\b\b\b%3.0f%%',ifile*100/nfiles);
end

% Supress Divide by Zero warnings in case of missing data points 
warning('off');

 % Return mean_val
mean_val = epr_struct(epr_acc./epr_count,min(start_time),max(stop_time),'Mean');
if min(min(epr_count)) == 0
  fprintf(1,'\nmean_val contains %.0f%% empty values (NaN).',sum(sum(isnan(mean_val.data)))*100/(180*288));
end

% Return mean_lat
if nargout > 1
  mean_lat = sum(epr_acc')'./sum(epr_count')';
  if min(sum(epr_count')) == 0
    fprintf(1,'\nmean_lat contains %.0f%% empty values (NaN).',sum(isnan(mean_lat))*100/180);
  end
end

fprintf('\n');