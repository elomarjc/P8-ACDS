function [ q ] = A2q( A )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

persistent q_old

a = A;
q_test = zeros(4,1);
q = zeros(4,1);
% Find i for largest denominater
q_test(4) = 0.5*sqrt(1+a(1,1)+a(2,2)+a(3,3));
q_test(1) = 0.5*sqrt(1+a(1,1)-a(2,2)-a(3,3));
q_test(2) = 0.5*sqrt(1-a(1,1)+a(2,2)-a(3,3));
q_test(3) = 0.5*sqrt(1-a(1,1)-a(2,2)+a(3,3));

dmax = [q_test(4), q_test(1), q_test(2), q_test(3)];

[nymax, i] = max(dmax);

if (i==1)
    q(4) = 0.5*sqrt(1+a(1,1)+a(2,2)+a(3,3));
    q(1) = 0.25*(a(2,3)-a(3,2)) / q(4);
    q(2) = 0.25*(a(3,1)-a(1,3)) / q(4);
    q(3) = 0.25*(a(1,2)-a(2,1)) / q(4);
elseif (i==2)
    q(1) = 0.5*sqrt(1+a(1,1)-a(2,2)-a(3,3));
    q(2) = 0.25*(a(1,2)+a(2,1)) / q(1);
    q(3) = 0.25*(a(1,3)+a(3,1)) / q(1);
    q(4) = 0.25*(a(2,3)-a(3,2)) / q(1); % Here Sidi has a sign error
elseif (i==3)
    q(2) = 0.5*sqrt(1-a(1,1)+a(2,2)-a(3,3));
    q(1) = 0.25*(a(1,2)+a(2,1)) / q(2);
    q(3) = 0.25*(a(2,3)+a(3,2)) / q(2);
    q(4) = 0.25*(a(3,1)-a(1,3)) / q(2);
elseif (i==4)
    q(3) = 0.5*sqrt(1-a(1,1)-a(2,2)+a(3,3));
    q(1) = 0.25*(a(1,3)+a(3,1)) / q(3);
    q(2) = 0.25*(a(2,3)+a(3,2)) / q(3);
    q(4) = 0.25*(a(1,2)-a(2,1)) / q(3);
end

% Return q
q = qunit(q);

% Remember the old quaternion
if (isempty(q_old))
    q_old = q;
elseif (q_old'*q < 0)
    q = -q;
end
q_old = q;
end
