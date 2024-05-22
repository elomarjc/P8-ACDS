% Calculates the flux comming from the earth reflection on
% the basis of the FOV meaned reflectivity matrix
% and the unit sun vector from earth to sun in ECEF
% and the zenith vector in ECEF.
%
% Parameters:
%    rs: Is a unit vector to the sun in ECEF.
%    zenith: Is a vector to the satellite in ECEF.
%    mrm: The FOV meaned reflectivity matrix.
%
% Output: The flux reflected to the satellite calculated 
%         from the FOV meaned refletivity matrix
%
% DOES NOT WORK!!
function flux = flux_ref(rs, zenith, mrm)

sun_flux = 1366.9;

%The projection of the sun flux vector onto the unit zenith
flux = sun_flux * dot(rs, zenith/norm(zenith));