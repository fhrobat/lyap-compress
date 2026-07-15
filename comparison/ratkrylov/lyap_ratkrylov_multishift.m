function [U, Y, nmatvecs, k] = lyap_ratkrylov_multishift(A,c,normc,options) 
    % lyap_ratkrylov_adaptive: Rational Krylov subspace method with multishift cg for the
    % Lyapunov equation  A*X + X*A' = c*c'
    %
    %
    %
    % Inputs:
    %   A       : n x n Symmetric Positive Definite matrix
    %   c       : n x 1 right-hand side vector
    %   options : struct with fields
    %               .xi         - array of poles
    %               .tol_solver – inner iterative solver tolerance
    %               .maxit      – inner iterative solver max iterations
    %
    % Outputs:
    %   U       : n x m 
    %   Y       : m x m  solution of the projected Lyapunov equation
    %   k       : number of  poles used
    %   nmatvecs: total inner-solver iterations accumulated



xi = options.xi;
k = length(xi);

% U orthonormal basis of Q(A, c, xi)
nmatvecs = 0;
switch options.variant
    case 'multishiftcg'
        [U, nmatvecs] = multishiftCG(A, c/normc, xi, options.tol_solver, options.maxit);
        U = orth(U,eps);
        if isnumeric(A) == 1
            ProjMat = U' * A * U;
        else
            ProjMat = U' * A(U);
        end

    case 'multishiftcg_rem'
        [U, nmatvecs] = multishiftCG_rem(A, c/normc, xi, options.tol_solver, options.maxit);
        U = orth(U,eps);
        if isnumeric(A) == 1
            ProjMat = U' * A * U;
        else
            ProjMat = U' * A(U);
        end
end

% solve the projected equation
Y = lyap(ProjMat, -(U'*c) * (U'*c)'); 
Y = (Y + Y')/2;

end