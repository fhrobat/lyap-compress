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
if nargin < 6
    sz = size(v1,1);
    Q = zeros(sz, m);
    Q(:,1) = v2;
    beta = zeros(m,1);
    alpha = zeros(m,1);

    for i = 1:m
        if i == 1
            [Q(:, i+1), alpha(i), beta(i)] = short_recurrence_Lanczos(mult, v1, v2, b);
            v1 = v2;
            v2 = Q(:, i+1);

        elseif i < m
            [Q(:, i+1), alpha(i), beta(i)] = short_recurrence_Lanczos(mult, v1, v2, beta(i-1));
            v1 = v2;
            v2 = Q(:, i+1);
        else
            [v, alpha(i), beta(i)] = short_recurrence_Lanczos(mult, v1, v2, beta(i-1));
            v1 = v2;
            v2 = v;

        end

    end
else
    sz = size(v1,1);
    Q = zeros(sz, m);
    Q(:,1) = v2;

    for i = 1:m
        if i == 1
            Q(:, i+1) = short_recurrence_Lanczos(mult, v1, v2, b, beta(1), alpha(1));
            v1 = v2;
            v2 = Q(:, i+1);
        elseif i < m
            Q(:, i+1) = short_recurrence_Lanczos(mult, v1, v2, beta(i-1), beta(i), alpha(i));
            v1 = v2;
            v2 = Q(:, i+1);
        else
            v = short_recurrence_Lanczos(mult, v1, v2, beta(i-1), beta(i), alpha(i));
            v1 = v2;
            v2 = v;

        end

    end
end

end