%Fits a polynomial to the meaned albedo data as a function of latitude.
%
% Param:
%   p_order: The order of the polynomial to fit.
%
% Output:
%   Coefficient vector with length of the polynomial order + 1
%   for use with polyval
function p_coef = ob_alb_poly(p_order)

%loading data file
albedo_data = load('../../../space_environment_emulation/albedo/toolbox/epr_data/2003/ga030101-031231.mat');

% polynomial fitting of data
p_order = 2;
p_coef = polyfit(-89.5:1:89.5,mean(albedo_data.data')*100,p_order);

