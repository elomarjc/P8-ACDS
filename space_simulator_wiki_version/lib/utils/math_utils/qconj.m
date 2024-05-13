function [ q_conj ] = qconj( q )
% The complex conjugated Quaternion
q_conj = [-q(1) -q(2) -q(3) q(4)]';
end

