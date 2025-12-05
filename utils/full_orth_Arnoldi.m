function [Q, v1, v2, H] = full_orth_Arnoldi(mult,m,v1, H)
%
% Function to compute m vectors of the Arnoldi basis starting from the vector v1
%
% Inputs: 
%   mult        function that handles matrix-vector product
%   m           number of vectors to compute
%   v1          starting vector
%   H           upper Hessenberg matrix H (optional). If provided, the algorithm does not
% compute the orthogonalization coefficients but takes them directly from H 
%
% Outputs: 
%   Q            Arnoldi basis
%   v1           last vector of Q 
%   v2           next basis vector
%   H            upper Hessenberg matrix of size m x m

if nargin < 4
    sz = size(v1,1);
    Q = zeros(sz,m);
    Q(:,1) = v1/norm(v1);
    H = zeros(m+1,m);
    for i = 1:m
        v = mult(Q(:,i));
        for j = 1:i
            H(j,i) = Q(:,j)'*v;
            v = v - H(j,i)*Q(:,j);
        end
        H(i+1,i) = norm(v);
        if i < m
            Q(:,i+1) = v/H(i+1,i);
        else 
            w = v/H(i+1,i);
        end
    end
    v1 = Q(:,m);
    v2 = w;
else
    sz = size(v1,1);
    Q = zeros(sz,m);
    Q(:,1) = v1/norm(v1);
    for i = 1:m
        v = mult(Q(:,i));
        for j = 1:i
            v = v - H(j,i)*Q(:,j);
        end
        if i < m
            Q(:,i+1) = v/H(i+1,i);
        else
            w = v/H(i+1,i);
        end
    end
    v1 = Q(:,m);
    v2 = w;
end

