% Finds the corresponding reflectivity on the basis
% of the FOV meaned reflectivity matrix and 
% a zenith vector in ECEF.
%
% Parameters:
%    zenith: Is a vector to the satellite in ECEF.
%    mrm: The FOV meaned reflectivity matrix.
%
% Output: The flux reflected to the satellite found
%         from the FOV meaned refletivity matrix
function ref = zenith2ref(zenith, mrm)

[z(1) z(2) z(3)]= cart2sph(zenith(1),zenith(2),zenith(3));

lat = round(89.5 + 180/pi*(z(2)) +1);
long = round(((179.375+180/pi*(z(1)))/(360/288)) + 1);

ref = mrm(lat, long);
