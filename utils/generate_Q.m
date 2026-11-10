function [Q, v1, v2, alpha, beta] = generate_Q(mult, m, v1, v2, b, beta, alpha)
%
% Function to extend Krylov basis by m vectors using Lanczos algorithm
%
% Inputs:
%   mult            function that handles matrix-vector product
%   m               number of vectors to compute
%   v2              starting vector (last vector of Lanczos basis)
%   v1              second-last vector
%   b               last beta coefficient
%   beta, alpha     beta and alpha coefficients if already known (optional)
%
% Outputs:
%   Q               the m new basis vectors (the first one is v2)
%   v1, v2          last vector of Q and next Lanczos vector
%   alpha, beta     Lanczos coefficients (beta(m) couples the last vector of Q and the next vector)
%
known = (nargin >= 6);
sz = size(v1,1);
Q = zeros(sz, m);
Q(:,1) = v2;
if ~known
    beta = zeros(m,1);
    alpha = zeros(m,1);
end

for i = 1:m
    % beta coefficient coupling v1 and v2
    if i == 1
        b_prev = b;
    else
        b_prev = beta(i-1);
    end

    if known
        v = short_recurrence_Lanczos(mult, v1, v2, b_prev, beta(i), alpha(i));
    else
        [v, alpha(i), beta(i)] = short_recurrence_Lanczos(mult, v1, v2, b_prev);
    end
    v1 = v2;
    v2 = v;

    if i < m
        Q(:, i+1) = v;
    end
end

end
