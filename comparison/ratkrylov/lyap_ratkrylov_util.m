function [U, Y, nmatvecs, k] = lyap_ratkrylov_util(A,c,normc,options) 


xi = options.xi;
k = length(xi);

% set rational Krylov to compute real basis
options.isreal = 1;
if isfield(options, 'precond') == 0
        options.precond = '';
end


% U orthonormal basis of Q(A, c, xi)
nmatvecs = 0;
switch options.variant
    case 'backslash'
        U = rational_krylov(A, c/normc, xi, options);
    case 'pcg'
        [U, nmatvecs] = rational_krylov_pcg(A, c/normc, xi, options);
    case 'multishiftcg'
        [U, nmatvecs] = multishiftCG(A, c/normc, xi, options.tol_solver, options.max_iter);
        U = orth(U,eps);

    case 'multishiftcg_rem'
        [U, nmatvecs] = multishiftCG_rem(A, c/normc, xi, options.tol_solver, options.max_iter);
        U = orth(U,eps);
end

if isnumeric(A) == 1
    ProjMat = U' * A * U;
else
    ProjMat = U' * A(U);
end

% solve the projected equation
Y = lyap(ProjMat, -(U'*c) * (U'*c)'); 
Y = (Y + Y')/2;


end