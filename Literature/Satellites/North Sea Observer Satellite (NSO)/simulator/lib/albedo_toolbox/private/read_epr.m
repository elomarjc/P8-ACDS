% READ_EPR Reads EP/TOMS epr data from file into an array. The return value
% is a struct containing the data and the start and stop times in julian
% date format. The type field specifies what kind of data is in the data
% field, e.g. average or raw.
%
% epr = read_epr(filename)
%
% where filename is a text string specifying the file to read.
%
% $Id: read_epr.m,v 1.1 2006/03/17 12:13:57 mnkr02 Exp $

function epr = read_epr(eprfile)

%%%%%%%%%%%%%%%%%%
% Init
%%%%%%%%%%%%%%%%%%

if nargin < 1
    % No file specified
    error('No file specified.');
end

[fid msg] = fopen(eprfile,'r');

if fid == -1
  error(msg);
end

%%%%%%%%%%%%%%%%%%
% Read header
%%%%%%%%%%%%%%%%%%

% Read date line
start_time = get_date(fid);
stop_time = start_time+1;

% Discard header line 2 and 3
junk = fgetl(fid);
junk = fgetl(fid);

%%%%%%%%%%%%%%%%%%
% Read data
%%%%%%%%%%%%%%%%%%

data = zeros(180,288);
fprintf(1,'%3.0f%%',0);

% Latitude from 0 to 179 deg in 1 deg steps (-89.5 to 89.5)
for ilat = 1:180
    % Read 11 lines with 25 values
    for line = 1:11
        % Read empty space
        junk = fscanf(fid,'%c',1);        
        % Read 25 values in 3-digit codes EMM = MM*10^E
        for value = 1:25
            data(ilat,(line-1)*25+value) = get_value(fid);    
        end
        % Read rest of line
        junk = fgetl(fid);
    end
    % Read 1 line with 13 values
    % Read empty space
    junk = fscanf(fid,'%c',1);        
    % Read 13 values in 3-digit codes EMM = MM*10^E (M.M*10 for epr) 
    for value = 1:13
        data(ilat,275+value) = get_value(fid);    
    end
    % Read rest of line
    junk = fgetl(fid);
    % Print status
    fprintf(1,'\b\b\b\b%3.0f%%',ilat*100/180);
end

fclose(fid);

epr = epr_struct(data,start_time,stop_time,'Raw');

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read 3-digit code and convert to double
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function value = get_value(fid)

e = fscanf(fid,'%c',1);
if e == ' '
    e = '0';
end
m = fscanf(fid,'%c',2);

if (e == '9' & m == '99')
    value = NaN;
else
	  % *0.1 removed for epr files. 0.01  added for fraction instead of
	  % percent
    value = str2double(strcat(m,'e',e))*0.01; 
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract Julian Date from header line 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function jd = get_date(fid);

% Read Day number
junk = fscanf(fid,'%c',10);

% Read Month
switch fscanf(fid,'%c',3);
  case 'Jan'
    M = 1;
  case 'Feb'
    M = 2;
  case 'Mar'
    M = 3;
  case 'Apr'
    M = 4;
  case 'May'
    M = 5;
  case 'Jun'
    M = 6;
  case 'Jul'
    M = 7;
  case 'Aug'
    M = 8;
  case 'Sep'
    M = 9;
  case 'Oct'
    M = 10;
  case 'Nov'
    M = 11;
  case 'Dec'
    M = 12;
  otherwise
    error('Error reading month in epr file');
end

junk = fscanf(fid,'%c',1);

% Read day
d1 = fscanf(fid,'%c',1);
if d1 == ' '
    d1 = '0';
end

d2 = fscanf(fid,'%c',1);

D = str2num([d1 d2]);

junk = fscanf(fid,'%c',2);

% Read year
Y = str2num(fscanf(fid,'%c',4));

% Read rest of line
junk = fgetl(fid);

% Calculate Julian Date (Cite: Wolfram Astronomy)
jd = d2jd(Y,M,D);

return

