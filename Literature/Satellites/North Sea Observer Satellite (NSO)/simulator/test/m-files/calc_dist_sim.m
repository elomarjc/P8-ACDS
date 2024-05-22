% When excecuting this file the mat file dist_sim.mat should be loaded
% first. This file calculates the mean, standard deviation and maximum of
% all the samples obtained from the monte carlo simulation in dist_sim.m

%%%%%% IMPORTANT LOAD dist_sim.mat %%%%%%%%
x=[];
for N=1:39
    for k=1:30
        x = [x sim_dist{N}{k}(:,:)'];
    end
end

mean_dist_x = mean(x(1,:));
mean_dist_y = mean(x(2,:));
mean_dist_z = mean(x(3,:));

mean_dist = [mean_dist_x mean_dist_y mean_dist_z];

std_dist_x = std(x(1,:));
std_dist_y = std(x(2,:));
std_dist_z = std(x(3,:));

std_dist = [std_dist_x std_dist_y std_dist_z];

max_dist_x = max(x(1,:));
max_dist_y = max(x(2,:));
max_dist_z = max(x(3,:));

max_dist = [max_dist_x max_dist_y max_dist_z];

min_dist_x = min(x(1,:));
min_dist_y = min(x(2,:));
min_dist_z = min(x(3,:));

min_dist = [min_dist_x min_dist_y min_dist_z];

save dist_mean_max_min_std.mat mean_dist std_dist max_dist min_dist; 
