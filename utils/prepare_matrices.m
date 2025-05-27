function [mult,c,options] = prepare_matrices(name, k)
%
% Function to prepare coefficient matrix, right-hand side and options (if
% needed)
%

options = [];

% prepare matrices for Example 1
if name == "4DLap"
    A = (k+1)^2*gallery('tridiag',k);
    A = kron(eye(k),A) + kron(A,eye(k));
    mult = @(v) A*v;
    c = zeros(k,1);
    for i = 1:k
        c(i) = 1/(0.5*sqrt(2*pi))*exp(-(i/(k+1)-1/2)^2/(2*0.5^2));
    end
    c = kron(c,c);
    options.eigmin = (k+1)^2 * (4 - 4 * cos(pi / (k+1)));
    options.eigmax = (k+1)^2 * (4 - 4 * cos(k * pi/ (k+1)));

% prepare matrices for Example 2
elseif name == "rail"
    switch k
        case 1
            load("Example2/rail_5177.mat");
        case 2
            load("Example2/rail_20209.mat");
        case 3
            load("Example2/rail_79841.mat");
    end
    p = dissect(E);
    R = chol(E(p,p));
    L = R';
    mult = @(v) -L\(A*((L')\v));
    c = - L\B(:,1);

% prepare matrices for Example 3
elseif name == "mor"
    switch k
        case 1
            load("Example3/mor_4813_10.mat");
        case 2
            load("Example3/mor_13551_10.mat");
        case 3
            load("Example3/mor_25872_10.mat");
        case 4
            load("Example3/mor_39527_10.mat");
    end
    L = diag(sqrt(diag(E)));
    A = - (L \ A) / L';
    A = (A + A')/2;
    mult = @(v) A * v;
    c = - L\B(:,1);

end


end