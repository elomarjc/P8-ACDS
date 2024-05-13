function [q] = dcmatrix_and_q(input)
% Direction Cosine Matrix computation and conversion to quaternion
% input=(O1;O2;O3) are the three unit vectors representing the axes of the 
% rotated coordinate system given in the reference coordinate system.

    O1=input(1:3);
    O2=input(4:6);
    O3=input(7:9);
    %O=[O1,O2,O3];
  
    % Make sure that the unit vectors are orthogonal
    %[Q,R]=qr(O);
    %O1=Q(1:3,1);
    %O2=Q(1:3,2);
    %O3=Q(1:3,3);
   
%     % Reference coordinate system axes:
%     I1=[1;0;0];
%     I2=[0;1;0]; 
%     I3=[0;0;1];

    % Direction Cosine Matrix
%     A = [O1'*I1 , O1'*I2 , O1'*I3;
%          O2'*I1 , O2'*I2 , O2'*I3;
%          O3'*I1 , O3'*I2 , O3'*I3];
    A = [O1 O2 O3];
    
    % Check if it belongs to SO(3) 
    %determinantA=det(A) % must be plus 1.
    %AAtrans = A*A'      % must be equal to ->
    %AtransA = A'*A      % this one.

    % Convert A to a quaternion
%     q = A_to_q(A);
      q = A2q(A);
    
    % Check if q is a valid rotation 
    %qnorm = q'*q       % must have norm 1
end