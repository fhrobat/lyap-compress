function [U, X, result] = lyap_compress(A,c,options) 
%
% Function to compute an approximation of the solution of the Lyapuov equation 
%               A * X + X * A = c * c'
% using the compress algorithm.
% The approximate solution is returned in a factorized form, i.e. U*X*U'

%Inputs: 
%   A                   real symmetric sparse matrix of size n x n or a function that handles matrix-vector product with A
%   c                   real vector of size n
%   options.tol         relative tolerance on the residual norm
%   options.xi          poles for the rational Krylov subspace (optional)
%   options.eigmin
%   options.eigmax      estimates of the extreme eigenvalues of A (optional)
%   options.maxmem      max number of vectors to keep in memory 
%   options.true_res    if set to 1, the true residual norm is computed and stored whener the estimate residual norm is computed
%
% Outputs:
%   U, X                matrices to form the approximate solution U*X*U'
%   result.iter         total number of Lanczos iterations
%   result.k            number of poles
%   result.m            set equal to maxmem - 2*k - 1 and indicates how often the compression is performed
%   result.errest       matrix with 2 columns showing the number of iterations (first column) and the relative estimate of  the residual norm (second column)
%   result.eigmin
%   result.eigmax       estimates of the extreme eigenvalues of A (or options.eigmin and options.eigmax if provided as inputs)
%   result.xi           list of poles

normc = norm(c);
flag = 0;

if ~isa(A,'function_handle')
		mult = @(v) A*v;
	else 
		mult = A;
end

if nargin < 3
		error("Please provide relative tolerance and maximum number of vectors to be kept in memory \n");

elseif isfield(options, 'tol') == 0
        error("Please provide relative tolerance \n")

elseif isfield(options, 'maxmem') == 0
        error("Please provide maximum number of vectors to be kept in memory \n")

elseif isfield(options, 'xi') == 0 && (isfield(options, 'eigmin') == 0 || isfield(options, 'eigmax') == 0)
    % if the poles are not provided and at leat one of the extreme eigevalues are not
    % provided, the first pass is performed employing full orthogonalization to estimate the eigenvalues
    flag = 1;
    if isfield(options, 'eigmin') == 1
        [Q, v1, v2, H] = full_orth_Arnoldi(mult, options.maxmem-1, c/normc);
        alpha = diag(H);
        beta = diag(H,-1);
        options.eigmin = eigs(spdiags([[beta(1:end-1);0],alpha,[0;beta(1:end-1)]], [-1,0,1], length(alpha), length(alpha)), 1, 'smallestabs');
    elseif isfield(options, 'eigmax') == 1
        [Q, v1, v2, H] = full_orth_Arnoldi(mult, options.maxmem-1, c/normc);
        alpha = diag(H);
        beta = diag(H,-1);
        options.eigmax = eigs(spdiags([[beta(1:end-1);0],alpha,[0;beta(1:end-1)]], [-1,0,1], length(alpha), length(alpha)), 1, 'largestabs');
    else
        [Q, v1, v2, H] = full_orth_Arnoldi(mult, options.maxmem-1, c/normc);
        alpha = diag(H);
        beta = diag(H,-1);
        options.eigmin = eigs(spdiags([[beta(1:end-1);0],alpha,[0;beta(1:end-1)]], [-1,0,1], length(alpha), length(alpha)), 1, 'smallestabs');
        options.eigmax = eigs(spdiags([[beta(1:end-1);0],alpha,[0;beta(1:end-1)]], [-1,0,1], length(alpha), length(alpha)), 1, 'largestabs');
    end
    options.eigmin = options.eigmin /10;
    options.eigmax = options.eigmax * 1.1;
    c0 = options.eigmax/options.eigmin;
    options.xi = poles_Zolotarev(options.eigmin, options.eigmax, options.tol/(c0 * 2));
elseif isfield(options, 'xi') == 0
    c0 = options.eigmax/options.eigmin;
    options.xi = poles_Zolotarev(options.eigmin, options.eigmax, options.tol/(c0 * 2));
end

if isfield(options, 'true_res') == 0
    options.true_res = 0;
end

k = length(options.xi);
if options.maxmem < 2*k
    error("The variable maxmem is too low to achieve the requested accuracy, please provide maxmem which is at least %d \n", 2*k+1)
elseif options.maxmem < (3*k-1) % if compression is too frequent, compress algorithm may suffer
    warning("Compression is performed every %d Lanczos iterations, please consider increasing maxmem \n", options.maxmem - 2*k-1)
else 
    fprintf("Compression is performed every %d Lanczos iterations \n", options.maxmem - 2*k -1)
end

options.m = options.maxmem - 2*k-1;

if flag == 0
    [U, X, iter, k, m, errest] = lyap_compress_util(mult,c, normc, options);
else
    [U, X, iter, k, m, errest] = lyap_compress_util(mult,c,normc, options, Q, v1, v2, alpha, beta);
end
result = [];
result.iter = iter;
result.k = k;
result.m = m;
result.errest = errest;
result.eigmin = options.eigmin;
result.eigmax = options.eigmax;
result.xi = options.xi;