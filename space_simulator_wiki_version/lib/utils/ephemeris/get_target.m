function [out] = get_target(q_e_i,R_sat,Re,min_elevation,target_list)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
km2m = 1e3;
load(target_list);
target_name = 0;
%% Constructing the target struct
n = size(targets,1);
for i = 1:n
    value_name{i} = targets{i,1};
    value_lat{i} = targets{i,2}(1);
    value_lon{i} = targets{i,2}(2);
    value_alt{i} = targets{i,2}(3);
end

name = 'name';
lat = 'lat';
lon = 'lon';
alt = 'alt';

s_targets = struct(name,value_name,lat,value_lat,lon,value_lon,alt,value_alt);

%% Find the target vector from the known elements and the satellite position
% Iterate through the target list
for k = 1:n
    r_target(:,k) = [cosd(s_targets(k).lat)*cosd(s_targets(k).lon)...
        cosd(s_targets(k).lat)*sind(s_targets(k).lon)...
        sind(s_targets(k).lat)]';
    
    R_target_e(:,k) = (s_targets(k).alt + km2m*Re) * r_target(:,k);
    R_target_i(:,k) = qRot(R_target_e(:,k),q_e_i);
    
    D_target(:,k) = R_sat - R_target_i(:,k);
    sat_inclination(k) = acosd((D_target(:,k)'* R_target_i(:,k))/(norm(D_target(:,k))*norm(R_target_i(:,k))));
end
inSight = find(abs(sat_inclination) < (90 - min_elevation));

if length(inSight) > 1
    for i = 1:length(inSight)
        elevation(i) = abs(sat_inclination(inSight(i)));
    end
    [M,I] = min(elevation);
    target_name = inSight(I);
    r_target = R_target_e(:,inSight(I));
elseif isempty(inSight)
    target_name = 0;
    r_target = [0 0 0]';
elseif ~isempty(inSight)
    target_name = inSight;
    r_target = R_target_e(:,inSight);
end
out = [target_name;r_target];
end