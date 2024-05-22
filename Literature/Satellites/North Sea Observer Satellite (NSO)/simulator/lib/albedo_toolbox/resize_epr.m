% RESIZE_EPR Resize epr struct data.
%
% newepr = resize_epr(oldepr,redfac)
%
% where redfacn is the reduction factor.
%
% $Id: resize_epr.m,v 1.1 2006/03/17 12:13:52 mnkr02 Exp $

function newepr = resize_epr(oldepr,redfac);

if redfac > 0
  newepr = epr_struct(resizem(oldepr.data,1/redfac),oldepr.start_time,oldepr.stop_time,oldepr.type);
else
  newepr = oldepr;
end

return