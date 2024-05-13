function [ p ] = qmult( q1,q2 )
% ======================================================================
%> @brief Quaternion multiplication
%>
%> @param q1: left
%> @param q2: right
%>
%> @retval p: result
% ======================================================================

assert(size(q1,1) == 4, 'qmult: argument 1 is of wrong dimension (not 4x1)');
assert(size(q2,1) == 4, 'qmult: argument 2 is of wrong dimension (not 4x1)');

p = [q1(4)*q2(1:3) + q2(4)*q1(1:3) + cross(q1(1:3),q2(1:3));
    q1(4)*q2(4) - dot(q1(1:3),q2(1:3))];
end

