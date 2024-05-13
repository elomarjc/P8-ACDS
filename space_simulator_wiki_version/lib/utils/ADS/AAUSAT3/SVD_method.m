function [output] = SVD_method(a1,b1_i,b1_s,a2,b2_i,b2_s,init)
tic_s=tic;
% The SVD_method solves Wahba's problem by Singular Value Decomposition,
% originally proposed by Markley.
%
%   Inputs:
%       a1 = Weight on first vector observation (usually 1/sigma^2)
%       a2 = Weight on second vector observation (usually 1/sigma^2)
%       b1_i = First vector observation in ECI.
%   	b1_s = First vector observation in SBRF.
%       b2_i = Second vector observation in ECI.
%   	b2_s = Second vector observation in SBRF.
%       init = Initial run indicator (1=true,0=false).
%
%   Output:
%       The quaternion based on the the optimal rotation matrix A.
%       Error indicator (1=quaternion cannot be trusted, 0=succes).
%
% Implementation by K. Vinther and K. F. Jensen, Aalborg University, 2010.

    persistent tictoc;

    B=a1*b1_s*(b1_i)'+a2*b2_s*(b2_i)';

    % Singular Value Decompositionn, B=U*S*V'
    [U,S,V]= svdmex(B);
    temp=det(U)*det(V);
    
    % Check covariance P (of the rotation angle error vector)
    s1=S(1,1);
    s2=S(2,2);
    s3=temp*S(3,3);
    P=U*diag([(s2+s3)^(-1),(s3+s1)^(-1),(s1+s2)^(-1)])*U';
    temp2=sum(P);
    P_sum=sum(temp2);
    
    % If the covariance is big then the attitude is unobservable
    if P_sum < 0.1
        % Calculation of optimal attitude matrix
        Aopt=U*diag([1,1,temp])*V';

        % Attitude matrix to quaternion conversion
        output(1:4)=A_to_q(Aopt);
        
        % invert quaternion
        %output(1:3)=-output(1:3);
    else
        output(1:4)=[0 0 0 1];
    end
    
    % output must be finite in order to pass to Matlab function blok
    if isinf(P_sum) == 1 || P_sum > 1000 || isnan(P_sum) == 1
        output(5)=1000;
    else
        output(5)=P_sum;
    end
    
    % Store algorithm execution time to workspace
    tictoc=[tictoc,toc(tic_s)];
    assignin('base','tictoc_svd',tictoc);
    if init==1
        tictoc=[];
    end
end