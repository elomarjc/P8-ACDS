% script for running diferent simulations

%logging
diary off; 
delete fov_simulation.log.txt
diary fov_simulation.log.txt
disp(['%%% logging output to fov_simulation.log.txt %%%'])

%simulation 
re = 6371010;

%with cos mask
for index = 0:4,
  save hmm_asdf index re
  clear all
  close all
  load hmm_asdf
  index
  re
  if index==4 ,
    or = re + 5e5 + 300
  else
    or = re + 5e5 + index * 2e5
  end
  force_recalc = true
  use_bogus=false
  force_reload = true
  display_figs = false
  use_flux_mask = true
  main_mean_alb_comp;
%    figure(3)
%     subplot(3,1,1);
%   image(a_toms*56);
%    subplot(3,1,2);
%   image(sum_matr*56);
%     subplot(3,1,3);
%   image(mrm*56);
  save(['with_cos_mask_' num2str(or)])
end

% no cos mask
for index=0:4,
  save hmm_asdf index re
  clear all
  close all
  load hmm_asdf
  if index==4,
    or = re + 5e5 + 300
  else
    or = re + 5e5 + index*2e5
  end
  force_recalc = true
  use_bogus=false
  force_reload = true
  display_figs = false
  use_flux_mask = false
  main_mean_alb_comp;
  save(['no_cos_mask_' num2str(or)])
end


diary off
