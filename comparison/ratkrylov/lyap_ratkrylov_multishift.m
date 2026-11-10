function [U, Y, nmatvecs, k] = lyap_ratkrylov_multishift(A,c,options)
    % lyap_ratkrylov_multishift: Rational Krylov subspace method with multishift CG for the
    % Lyapunov equation  A*X + X*A' = c*c'
    %
    % The basis of the rational Krylov subspace Q(A, c, xi) is computed by
    % solving all the shifted linear systems (A - xi_j I) x_j = c with multishift CG.
    %
    % Inputs:
    %   A       : n x n Symmetric Positive Definite matrix (or function handle for the product with A)
    %   c       : n x 1 right-hand side vector
    %   options : struct with fields
    %               .xi         - array of poles
    %               .variant    - 'multishiftcg' | 'multishiftcg_rem'
    %               .tol_solver – inner iterative solver tolerance
    %               .maxit      – inner iterative solver max iterations
    %
    % Outputs:
    %   U       : n x m orthonormal basis
    %   Y       : m x m  solution of the projected Lyapunov equation
    %   nmatvecs: total number of matrix-vector products with A
    %   k       : number of poles used

normc = norm(c);
xi = options.xi;
k = length(xi);

% U orthonormal basis of Q(A, c, xi)
nmatvecs = 0;
switch options.variant
    case 'multishiftcg'
        [U, nmatvecs] = multishiftCG(A, c/normc, xi, options.tol_solver, options.maxit);
    case 'multishiftcg_rem'
        [U, nmatvecs] = multishiftCG_rem(A, c/normc, xi, options.tol_solver, options.maxit);
    otherwise
        error("Unknown variant %s", options.variant);
end

% each multishift CG iteration requires one matrix-vector product with A;
% the iteration counter of multishiftCG(_rem) starts from 1 (multishiftCG_rem
% returns one counter per shift)
nmatvecs = max(nmatvecs) - 1;

% orthogonalize the computed vectors (tolerance equal to machine precision)
U = orth(U,eps);
if isnumeric(A) == 1
    ProjMat = U' * A * U;
else
    ProjMat = U' * A(U);
end
ProjMat = (ProjMat + ProjMat')/2;

% solve the projected equation
Y = lyap(ProjMat, -(U'*c) * (U'*c)');
Y = (Y + Y')/2;
end
