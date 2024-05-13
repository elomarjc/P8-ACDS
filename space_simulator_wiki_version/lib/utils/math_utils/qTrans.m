function [ T ] = qTrans( q )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
T = -0.5*[-q(4) -q(3)  q(2);
        q(3)  -q(4) -q(1);
        -q(2)  q(1)  -q(4)];
end

