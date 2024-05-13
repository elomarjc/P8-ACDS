function [q] = A_to_q(A)
% Attitude matrix to quaternion conversion
% http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixT
% oQuaternion/index.htm
%     Trace = A(1,1) + A(2,2) + A(3,3);
% 
%     if Trace > 0
%         S  = sqrt(Trace+1.0) * 2;
%         q4 = 0.25 * S;
%         q1 = (A(3,2) - A(2,3)) / S;
%         q2 = (A(1,3) - A(3,1)) / S;  
%         q3 = (A(2,1) - A(1,2)) / S;
%     elseif (A(1,1) > A(2,2)) && (A(1,1) > A(3,3))
%         S  = sqrt(1.0 + A(1,1) - A(2,2) - A(3,3)) * 2;
%         q4 = (A(3,2) - A(2,3)) / S; 
%         q1 = 0.25 * S;
%         q2 = (A(1,2) + A(2,1)) / S;  
%         q3 = (A(1,3) + A(3,1)) / S; 
%     elseif A(2,2) > A(3,3)
%         S  = sqrt(1.0 + A(2,2) - A(1,1) - A(3,3)) * 2;
%         q4 = (A(1,3) - A(3,1)) / S;
%         q1 = (A(1,2) + A(2,1)) / S;
%         q2 = 0.25 * S;  
%         q3 = (A(2,3) + A(3,2)) / S; 
%     else
%         S  = sqrt(1.0 + A(3,3) - A(1,1) - A(2,2)) * 2;
%         q4 = (A(2,1) - A(1,2)) / S;
%         q1 = (A(1,3) + A(3,1)) / S;
%         q2 = (A(2,3) + A(3,2)) / S;  
%         q3 = 0.25 * S; 
%     end
%     
%     q = [q1 q2 q3 q4]';

a = A;
q_test = zeros(4,1);
q = zeros(4,1);
% Find i for largest denominater
q_test(4) = 0.5*sqrt(1+a(1,1)+a(2,2)+a(3,3));
q_test(1) = 0.5*sqrt(1+a(1,1)-a(2,2)-a(3,3));
q_test(2) = 0.5*sqrt(1-a(1,1)+a(2,2)-a(3,3));
q_test(3) = 0.5*sqrt(1-a(1,1)-a(2,2)+a(3,3));
% dmax = [q(4), a(1,1), a(2,2), a(3,3)];
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

end