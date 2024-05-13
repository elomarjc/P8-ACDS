function [ q_inv ] = qinv( q )
% ======================================================================
%> @brief Quaternion inverse
%>
%> @param q: quaternion to be inverted
%>
%> @retval p_inv: result
% ======================================================================
assert(size(q,1) == 4, 'qinv: argument 1 is of wrong dimension (not 1x4)');

q_inv = qconj(q)/sqrt(dot(q,q));
end

