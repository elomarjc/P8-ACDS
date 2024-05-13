function q = qmethod(b, r, a, qold)
% The q-method (with some quaternion direction fixes)
% Inputs
% b - input vectors
% r - reference frame
% a - weights
% qold - old quaternion (to preserve direction)
%
% Output
% q - the quaternion which maps b’s into r’s
% n is equal to the number of columns of b
n=size(b,2);
% Initialization of B%
B=zeros(3,3);
% Calculate B = sum(a*b*r’)
for k=1:n %this loop is determined for wahba problem%
B = B+a(k) * ( (b(:,k)*((r(:,k))')) );%this equation can be found page 73 of the ADS report from last year%
end
% Calculate z
z = [   (B(2,3)-B(3,2));
        (B(3,1)-B(1,3));
        (B(1,2)-B(2,1))     ];
% Trace of b
tr = trace(B);
% Calculate k
K = [(B+B'-eye(3)*tr) z ; z' tr];

% Find the eigenvector
[vec,val] = eig(K, 'nobalance');
% Result is not at hand
[lambdamax,i] = max(diag(val));
q=vec(:,i);
% Due to the eigenvector calculation, the quaternion
% may be ‘‘flipped'', which we try to correct using
% the difference between now and the old quaternion
err=qdiff(qold,q);
if err >= pi/2;
% Difference greater than 2pi, flip
q=-q;
end
%The qdiff function is used to find the difference between two quaternion, which is used for
%correcting the direction.
function angle = qdiff(q1, q2)
q2(1:4) = [ -q2(1:3);
q2(4) ];
A= [ q1(4) q1(3) -q1(2) q1(1);
-q1(3) q1(4) q1(1) q1(2);
q1(2) -q1(1) q1(4) q1(3);
-q1(1) -q1(2) -q1(3) q1(4) ];
q = A*q2;
angle = 2*acos(q(4));