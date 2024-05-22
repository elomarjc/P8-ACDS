% REPLACE_NAN Replaces NaN values with annual mean or specified
% reflectivity data.
%
% new_epr = replace_nan(main_epr, param)
%
% main_epr is the reflectivity data in which undefined values will be
% replaced. If param is a an epr struct, data from param is copied to
% main_epr. If param is a string it is used as a path to the library of
% annual reflectivity data. The libaray should contain a directiry for each
% year of data, and the files should must be of the structure'gaYYMMDD.mat'
% and the annual mean data should be 'gaYY0101-YY1231.mat'. Per defalt,
% i.e. only main_epr is specified, the path is searched for the annual
% mean.
%
% $Id: replace_nan.m,v 1.1 2006/03/17 12:13:52 mnkr02 Exp $

function new_epr = replace_nan(main_epr, param);

persistent nan_epr;

% Check parameters
if nargin > 1 && isstruct(param)
    index = find(isnan(main_epr.data));
    main_epr.data(index) = param.data(index);
    msg = 'parameter.';
else
  year = jd2d(main_epr.start_time);
  if ~(isfield(nan_epr,'data') && year == jd2d(nan_epr.start_time))
    if year > 2000
      fyear = year - 2000;
    else
      fyear = year - 1900;
    end
    filename = sprintf('ga%02.0f0101-%02.0f1231.mat',fyear,fyear);
    file = which(filename);
    if isempty(file) && exist('param')
      filename = fullfile(param,num2str(year),filename);
      msg = filename;
    else
      msg = file;
    end
    nan_epr = load(filename);
  else
    msg = 'pre-loaded data.';
  end
  index = find(isnan(main_epr.data));
  main_epr.data(index) = nan_epr.data(index);
end

new_epr = main_epr;
if size(index,1) > 0
  disp(['replace_nan.m: ' num2str(size(index,1)) ' NaN replacements from ' msg]);
  new_epr.type = 'Mean supported raw';
end
