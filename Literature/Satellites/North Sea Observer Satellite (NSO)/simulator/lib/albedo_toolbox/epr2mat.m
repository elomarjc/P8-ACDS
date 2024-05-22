% EPR2MAT Convert TOMS reflectivity data to MAT file.
%
% epr2mat(filename,outdir)
%
% Filename is a string of either a single file or wildcard name. Output
% files are in outdir or same path as input if not specified. Only '.epr'
% files are read. Data should be loaded with 'epr = load(filename)', which
% will load data into an epr struct.
%
% $Id: epr2mat.m,v 1.1 2006/03/17 12:13:51 mnkr02 Exp $

function epr2mat(file, outdir);

if nargin < 2
  if isdir(file)
    outdir = file;
  else
    outdir = fileparts(file);
  end
end

files = dir(file);
n = size(files,1);

if n == 0
  error('File not found.');
end

for i=1:n
  fprintf('\rFile %i of %i: ',i,n);
  [pathstr,name,ext,versn] = fileparts(files(i).name);
  if strcmp(ext,'.epr')
    epr = read_epr(files(i).name);
    data = epr.data;
    start_time = epr.start_time;
    stop_time = epr.stop_time;
    type = epr.type;
    [pathstr,name,ext,versn] = fileparts(files(i).name);
    filename = fullfile(outdir, [name '.mat']);
    save(filename,'data','start_time','stop_time','type');
  end
end

fprintf('\n');
