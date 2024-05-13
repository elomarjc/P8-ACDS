% This file calculates the mean and standard deviation from the 
% monte-carlo simulation of the angular velocity controller. 
% Either the ang_test_mont.m need to be run first or the file
% ang_mont_estimation.mat (or ang_mont_no_estimation.mat) must
% be loaded as this file assumes the variable y_mont exists.

N_stat = size(y_mont,2);

ang_vel_all = [];

for k_stat=1:N_stat
  ang_vel_all = [ang_vel_all y_mont{k_stat}(:,1:3)'];
end

ang_mean_stat = mean(ang_vel_all,2)
ang_std_stat = std(ang_vel_all,0,2)
ang_max_stat = max(ang_vel_all,[],2)

