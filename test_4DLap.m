% First experiment (4D Laplacian) in the paper A. Casulli, F. Hrobat, D. Kressner "Lanczos with
% compression for symmetric Lyapunov equations", testing algorithms for 
% the case of 4D Laplacian


% results will contain cells with matrix dimension, compress results, two-pass results and restart results 
results = []; 

% j is the square root of N

for j= [floor(10^2*sqrt(18)),floor(10^2 * sqrt(36)), floor(10^2 * sqrt(72)), 1200] 

addpath(genpath('./'))

% generate function handle for mult, right-hand side and eigenvalues
name = "4DLap";
[mult,c,options] = prepare_matrices(name, j); 

% normalize by norm(c)
normc = norm(c);
c = c/normc;
mult = @(v) mult(v) / (normc)^2;
options.eigmin = options.eigmin / (normc^2);
options.eigmax = options.eigmax / (normc^2);
disp(name)
disp(j)

% set options
tol = 1e-6;
options.maxmem = 120;
options.tol = tol;
options_true_res = 0;

% compress
tic
[U, X, result_compress] = lyap_compress(mult,c,options);
result_compress.time = toc;

% true residual norm for compress
[~, Rsx] = qr([mult(U) * X,  U, c], 'econ');
[~, Rdx] = qr([U, mult(U) * X, -c], 'econ');
result_compress.true_res = norm(Rsx * Rdx', 'fro')/norm(c)^2;

% two-pass
tic
[U, X, result_two_pass] = two_pass_lanc(mult,c,options);
result_two_pass.time = toc;

% true residual norm for two-pass
[~, Rsx] = qr([mult(U) * X,  U, c], 'econ');
[~, Rdx] = qr([U, mult(U) * X, -c], 'econ');
result_two_pass.true_res = norm(Rsx * Rdx', 'fro')/norm(c)^2;

% restart
A = [];
A.multiply = @(v) -mult(v);

param = [];
param.max_restarts = 11000/ options.maxmem;
param.memory_max = options.maxmem; % set same memory costraint
param.norm = 'fro';
param.tol_res = tol;
param.tol_comp = tol*1e-4; % default tolerance
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
results{end+1} = {j^2, result_compress, result_two_pass, result_rest};


end
