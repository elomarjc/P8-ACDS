%This file runs simulations of the disturbance torques the rotation of the
%satellite is 90 degrees around the z-axis. The file runs 10 simulations
%one for every 10 degrees and the workspace variables are saved in mat
%files.

clear all;
clc;

theta=deg2rad(0);
sim('uncertainty_disturbance_torques')
save udt_0.mat

clear all;
clc;

theta=deg2rad(10);
sim('uncertainty_disturbance_torques')
save udt_10.mat

clear all;
clc;

theta=deg2rad(20);
sim('uncertainty_disturbance_torques')
save udt_20.mat

clear all;
clc;

theta=deg2rad(30);
sim('uncertainty_disturbance_torques')
save udt_30.mat

clear all;
clc;

theta=deg2rad(40);
sim('uncertainty_disturbance_torques')
save udt_40.mat

clear all;
clc;

theta=deg2rad(50);
sim('uncertainty_disturbance_torques')
save udt_50.mat

clear all;
clc;

theta=deg2rad(60);
sim('uncertainty_disturbance_torques')
save udt_60.mat

clear all;
clc;

theta=deg2rad(70);
sim('uncertainty_disturbance_torques')
save udt_70.mat

clear all;
clc;

theta=deg2rad(80);
sim('uncertainty_disturbance_torques')
save udt_80.mat

clear all;
clc;

theta=deg2rad(90);
sim('uncertainty_disturbance_torques')
save udt_90.mat