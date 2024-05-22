% MASK Make masked array from two equal arrays.
%
% res = mask(epr_data,mask,contrast)
%
% epr_data is the reflectivity "background" data. mask is an array of
% logical values of visible indices. contrast [0..1] specifies the contrast of
% visible and not visible data points in epr_data.
%
% $Id: mask.m,v 1.1 2006/03/17 12:13:57 mnkr02 Exp $

function res = mask(epr_data,mask,contrast);

if nargin < 3
	contrast = 0.3;
end

res = (~mask*contrast+mask).*epr_data;