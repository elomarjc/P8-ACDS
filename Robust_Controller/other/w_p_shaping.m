close all
clear all
clc
s = tf('s');
% This file aims to translate the different plant requirements into wp and wt, such that
% robust contorl theory can be applied

%% Example 1: Steady state error
% Plant with 50% ss err:

L = 10/(s+10);

bode(L)