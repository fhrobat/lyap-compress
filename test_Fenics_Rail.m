% Second experiment (Model order reduction: example 1) in the paper A. Casulli, F. Hrobat, D. Kressner "Lanczos with
% compression for symmetric Lyapunov equations", testing algorithms for the
% Model Order Reduction problem FEniCS Rail


% results will contain cells with matrix dimension, compress results, two-pass results and restart results
results = [];

for j=1:3

addpath(genpath('./'))

% generate function handle for mult and right-hand side (options is empty)
name = "rail";
[mult,M,c,options] = prepare_matrices(name, j);
L = M{1};
F = M{2}; %the coefficient matrix is equal to -L\(F/(L'))
disp(name)
disp(j)

% set options
tol = 1e-3;
options.maxmem = 120;
options.tol = tol;
options.true_res = 0;

%%% compress
tic
[U, X, result_compress] = lyap_compress(mult,c,options);
result_compress.time = toc;

% true residual norm for compress
[~, Rsx] = qr([mult(U) * X,  U, c], 'econ');
[~, Rdx] = qr([U, mult(U) * X, -c], 'econ');
result_compress.true_res = norm(Rsx * Rdx', 'fro')/norm(c)^2;
%%%

%%% two-pass
tic
[U, X, result_two_pass] = two_pass_lanc(mult,c,options);
result_two_pass.time = toc;

% true residual norm for two-pass
[~, Rsx] = qr([mult(U) * X,  U, c], 'econ');
[~, Rdx] = qr([U, mult(U) * X, -c], 'econ');
result_two_pass.true_res = norm(Rsx * Rdx', 'fro')/norm(c)^2;
%%%

%%% restart
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
%%%

%%%%%extended
param = [];
param.max_restarts = 11000/ options.maxmem;
param.memory_max = options.maxmem; % set same memory costraint
param.norm = 'fro';
param.tol_res = tol;
param.tol_comp = tol*1e-4; % default tolerance
param.verbose = 0;


ek_param = [];
ek_param.solver = 'bcg';                                        % options include bcg, bfom, and bgmres
ek_param.bcg_variant = 'bcgrq';
ek_param.symmetric = true;                                      % specification for struct_make
ek_param.max_iterations = 1000;
ek_param.max_restarts = param.max_restarts;
ek_param.verbose = 0;
ek_param.tol_res = tol*1e-2;

A = [];
A.multiply =@(x) -mult(x);
Astruct = struct( 'solve', @(x) L' * block_solve(F, L*x, ek_param),...
            'multiply', @(x) A.multiply(x),...
            'solver', ek_param.solver);

%%% cg without preconditioning
tic
param = extended_ksm_mat_eq(Astruct,[],c,[],param);
param.time = toc;
if ~isempty(param.S)
    [~, Rsx] = qr([A.multiply(param.Xleft)*param.S,  param.Xleft, c], 'econ');
    [~, Rdx] = qr([param.Xleft, A.multiply(param.Xleft)*(param.S)', c], 'econ');
    param.true_res = norm(Rsx * Rdx', 'fro')/norm(c)^2;
else
    [~, Rsx] = qr([A.multiply(param.Xleft),  param.Xleft, c], 'econ');
    [~, Rdx] = qr([param.Xleft, A.multiply(param.Xleft), c], 'econ');
    param.true_res = norm(Rsx * Rdx', 'fro')/norm(c)^2;
end

result_extended_plaincg = param;
%%%

param = [];
param.max_restarts = 11000/ options.maxmem;
param.memory_max = options.maxmem; % set same memory costraint
param.norm = 'fro';
param.tol_res = tol;
param.tol_comp = tol*1e-4; % default tolerance
param.verbose = 0;


%%%cg with ichol
tic

LL = ichol(-F);
ek_param.precond_left = -(LL);
ek_param.precond_right = -(LL)';
Astruct = struct( 'solve', @(x) L' * block_solve(F, L*x, ek_param),...
            'multiply', @(x) A.multiply(x),...
            'solver', ek_param.solver);

param = extended_ksm_mat_eq(Astruct,[],c,[],param);
param.time = toc;
if ~isempty(param.S)
    [~, Rsx] = qr([A.multiply(param.Xleft)*param.S,  param.Xleft, c], 'econ');
    [~, Rdx] = qr([param.Xleft, A.multiply(param.Xleft)*(param.S)', c], 'econ');
    param.true_res = norm(Rsx * Rdx', 'fro')/norm(c)^2;
else
    [~, Rsx] = qr([A.multiply(param.Xleft),  param.Xleft, c], 'econ');
    [~, Rdx] = qr([param.Xleft, A.multiply(param.Xleft), c], 'econ');
    param.true_res = norm(Rsx * Rdx', 'fro')/norm(c)^2;
end

result_extended_icholcg = param;
%%%


param = [];
param.max_restarts = 11000/ options.maxmem;
param.memory_max = options.maxmem; % set same memory costraint
param.norm = 'fro';
param.tol_res = tol;
param.tol_comp = tol*1e-4; % default tolerance
param.verbose = 0;


%%%precomputed chol
tic
LL = chol(-F);
LL = LL';
A.solve =@(x) -L'*(LL'\(LL\(L*x)));
A.solver = 'backslash';

param = extended_ksm_mat_eq(A,[],c,[],param);
param.time = toc;
if ~isempty(param.S)
    [~, Rsx] = qr([A.multiply(param.Xleft)*param.S,  param.Xleft, c], 'econ');
    [~, Rdx] = qr([param.Xleft, A.multiply(param.Xleft)*(param.S)', c], 'econ');
    param.true_res = norm(Rsx * Rdx', 'fro')/norm(c)^2;
else
    [~, Rsx] = qr([A.multiply(param.Xleft),  param.Xleft, c], 'econ');
    [~, Rdx] = qr([param.Xleft, A.multiply(param.Xleft), c], 'econ');
    param.true_res = norm(Rsx * Rdx', 'fro')/norm(c)^2;
end

result_extended_backslash = param;
%%%

%%%%%


% append results
res.name = "Fenics_Rail";
res.size = size(M,1);
res.compress =  result_compress;
res.two_pass =  result_two_pass;
res.extended_plaincg =  result_extended_plaincg;
res.extended_icholcg =  result_extended_icholcg;
res.extended_backslash =  result_extended_backslash;
res.rest =  result_rest;
results{end+1} = res;
end
