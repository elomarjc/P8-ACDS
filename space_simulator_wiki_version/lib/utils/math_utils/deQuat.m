function [ u,theta ] = deQuat( q )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

q_hat = q(1:3);
q_real = q(4);

theta = 2*acos(q_real);
v = (q_hat/sin(acos(q_real)));
u = v/norm(v);
end

