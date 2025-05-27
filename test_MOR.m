% Third experiment (Model order reduction: example 2) in the paper A. Casulli, F. Hrobat, D. Kressner "Lanczos with
% compression for symmetric Lyapunov equations", testing algorithms for the
% Model Order Reduction problem taken from Experiment 3 in P. Benner, D. Palitta, J. Saak 
% "On an integrated Krylov‑ADI solver for large‑scale Lyapunov equations"


% results will contain cells with matrix dimension, compress results, two-pass results and restart results
results = [];

for j=1:4

addpath(genpath('./'))

% generate function handle mult and right-hand side (options is empty)
name = "mor";
[mult,c,options] = prepare_matrices(name, j);
disp(name)
disp(j)

% set options
tol = 1e-3;
options.maxmem = 120;
options.tol = tol;
options.true_res = 0;

% compress
tic
[U, X, result_compress] = lyap_compress(mult,c,options);
result_compress.time = toc;

% true residual norm for compress
[~, Rsx] = qr([mult(U) * X,  U, c], 'econ');
[~, Rdx] = qr([U, mult(U) * X, -c], 'econ');
result_compress.true_res = norm(Rsx * Rdx', 'fro')/norm(c)^2

% two-pass
tic
[U, X, result_two_pass] = two_pass_lanc(mult,c,options);
result_two_pass.time = toc;

% true residual norm for two-pass
[~, Rsx] = qr([mult(U) * X,  U, c], 'econ');
[~, Rdx] = qr([U, mult(U) * X, -c], 'econ');
result_two_pass.true_res = norm(Rsx * Rdx', 'fro')/norm(c)^2

% restart
A = [];
A.multiply = @(v) -mult(v);


param = [];
param.max_restarts = 11000/ options.maxmem;
param.memory_max = options.maxmem; % set same memory costraint
param.norm = 'fro';
param.tol_res = tol;
param.tol_comp = eps; % tolerance set to eps
param.verbose = 0;
tic
param = restarted_ksm_mat_eq(A,[],c,[],param);
param.time_rest = toc;

% true residual norm for restart
if ~isempty(param.S)
    [~, Rsx] = qr([A.multiply(param.Xleft)*param.S,  param.Xleft, c], 'econ');
    [~, Rdx] = qr([param.Xleft, A.multiply(param.Xleft)*(param.S)', c], 'econ');
    param.true_res = norm(Rsx * Rdx', 'fro')/norm(c)^2;
else
    [~, Rsx] = qr([A.multiply(param.Xleft),  param.Xleft, c], 'econ');
    [~, Rdx] = qr([param.Xleft, A.multiply(param.Xleft), c], 'econ');
    param.true_res = norm(Rsx * Rdx', 'fro')/norm(c)^2;
end

% result for restart
result_rest = param;

% append results
results{end+1} = {size(c,1), result_compress, result_two_pass, result_rest};

end
