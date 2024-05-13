function [ v ] = qRot( u,q )
% ======================================================================
%> @brief Rotation of vector by quaternion
%>
%> @param u: vector to be rotated [3x1]
%> @param q: roatation quaternion
%>
%> @retval v: result
% ======================================================================
assert(size(u,1) == 3, 'qrot: argument 1 is of wrong dimension (not 3x1)');
assert(size(q,1) == 4, 'qmult: argument 2 is of wrong dimension (not 4x1)');

v_hat = qmult(qmult(qinv(q),[u;0]),q);
v = v_hat(1:3);

assert(size(v,1) == 3, 'v: wrong output dimension (not 3x1)');
end

