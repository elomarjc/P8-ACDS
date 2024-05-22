% EPR_STRUCT Create EPR struct from parameters. The struct contains the following
% fields:
%    - data (reflectivity data)
%    - start_time (Julian Date of start time)
%    - stop_time (Julian Date of stop time)
%    - type (String specifying data description, e.g. Raw or Mean)
%
% epr = epr_struct(data,start_time,stop_time,type)
%
% $Id: epr_struct.m,v 1.1 2006/03/17 12:13:51 mnkr02 Exp $

function epr = epr_struct(data,start_time,stop_time,type);

epr.data = data;
epr.start_time = start_time;
epr.stop_time = stop_time;
epr.type = type;
