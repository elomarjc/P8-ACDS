% ALBEDO_WRAPPER Simulink wrapper for ALBEDO.
%
% a = albedo_wrapper(sat,sun,param,redfac,eprlib);
%
% sat and sun are the vectors to the Earth and Sun in Earth Centered Earth
% Fixed coordinates, respectively. Param is either an epr struct or a
% julian date. redfac is the reduction factor, which is used to reduce the
% resolution of the epr data. eprlib is the path to the reflectivity data,
% only needed if param is a Julian Date.
%
% $Id: albedo_onboard_wrapper.m,v 1.1 2006/03/17 12:13:50 mnkr02 Exp $

function a = albedo_onboard_wrapper(sat,sun,param,redfac,eprlib);

CONST.AM0 = 1366.9;
persistent epr;

if isstruct(param)
  if ~isstruct(epr) || epr.start_time ~= param.start_time || epr.stop_time ~= param.stop_time
    epr = resize_epr(param,redfac);
  end
    a = albedo(sat,sun,epr);
else  
  if ~(isfield(epr,'data') && param >= epr.start_time && param <= epr.stop_time)
    epr = load(jd2file(param,eprlib));
    epr = replace_nan(epr,eprlib);
    epr = resize_epr(epr,redfac);
    % Check if data is unavailable. If so, earliest data is selected, and
    % is valid from current time to stop time.
    if param < epr.start_time
      epr.start_time = param;
    end
  end
  a = albedo(sat,sun,epr);
end
