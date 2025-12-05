% First experiment (4D Laplacian) in the paper A. Casulli, F. Hrobat, D. Kressner "Lanczos with
% compression for symmetric Lyapunov equations", testing algorithms for 
% the case of 4D Laplacian


% results will contain cells with matrix dimension and results for each method 
results = []; 

% j is the square root of N

for j= [floor(10^2*sqrt(18)),floor(10^2 * sqrt(36)), floor(10^2 * sqrt(72)), 1200] 

addpath(genpath('./'))

% generate function handle for mult, right-hand side and eigenvalues
name = "4DLap";
[mult,A,c,options] = prepare_matrices(name, j); 

% normalize by norm(c)
normc = norm(c);
c = c/normc;
mult = @(v) mult(v) / (normc)^2;
options.eigmin = options.eigmin / (normc^2);
options.eigmax = options.eigmax / (normc^2);
A = A / (normc^2);
disp(name)
disp(j)

% set options
tol = 1e-6;
options.maxmem = 120;
options.tol = tol;
options.true_res = 0;
options.isreal = 1;

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


%%%ratkrylov backslash
options.variant = 'backslash';
tic
[U, X, result_ratkrylov_backslash] =  lyap_ratkrylov(A,c,options);
result_ratkrylov_backslash.time = toc;

% true residual norm for ratkrylov backslash
[~, Rsx] = qr([mult(U) * X,  U, c], 'econ');
[~, Rdx] = qr([U, mult(U) * X, -c], 'econ');
result_ratkrylov_backslash.true_res = norm(Rsx * Rdx', 'fro')/norm(c)^2;
%%%


%%%ratkrylov plain cg
options.variant = 'pcg';
options.tol_solver = 1e-6;
tic
[U, X, result_ratkrylov_plaincg] =  lyap_ratkrylov(A,c,options);
result_ratkrylov_plaincg.time = toc;

% true residual norm for ratkrylov plain cg
[~, Rsx] = qr([mult(U) * X,  U, c], 'econ');
[~, Rdx] = qr([U, mult(U) * X, -c], 'econ');
result_ratkrylov_plaincg.true_res = norm(Rsx * Rdx', 'fro')/norm(c)^2;
%%%



%%%ratkrylov ichol pcg
options.variant = 'pcg';
options.precond = 'ichol';
options.tol_solver = 1e-6;
tic
[U, X, result_ratkrylov_icholcg] =  lyap_ratkrylov(A,c,options);
result_ratkrylov_icholcg.time = toc;

% true residual norm for ratkrylov ichol pcg
[~, Rsx] = qr([mult(U) * X,  U, c], 'econ');
[~, Rdx] = qr([U, mult(U) * X, -c], 'econ');
result_ratkrylov_icholcg.true_res = norm(Rsx * Rdx', 'fro')/norm(c)^2;
%%%
options = rmfield(options, 'tol_solver');


%%%ratkrylov mutlishiftcg
options.max_iter = 10000;
options.variant = 'multishiftcg';
tic
[U, X, result_ratkrylov_multishiftcg] =  lyap_ratkrylov(A,c,options);
result_ratkrylov_multishiftcg.time = toc;

% true residual norm for ratkrylov mutltishiftcg
[~, Rsx] = qr([mult(U) * X,  U, c], 'econ');
[~, Rdx] = qr([U, mult(U) * X, -c], 'econ');
result_ratkrylov_multishiftcg.true_res = norm(Rsx * Rdx', 'fro')/norm(c)^2;
%%%

%%%ratkrylov mutlishiftcg with removal
options.variant = 'multishiftcg_rem';
tic
[U, X, result_ratkrylov_multishiftcg_rem] =  lyap_ratkrylov(A,c,options);
result_ratkrylov_multishiftcg_rem.time = toc;

% true residual norm for ratkrylov mutlishiftcg with removal
[~, Rsx] = qr([mult(U) * X,  U, c], 'econ');
[~, Rdx] = qr([U, mult(U) * X, -c], 'econ');
result_ratkrylov_multishiftcg_rem.true_res = norm(Rsx * Rdx', 'fro')/norm(c)^2;
%%%


%%%%%extended
A = -A;
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


Astruct = struct_make(A, ek_param);

%%% cg without preconditioning
tic
param = extended_ksm_mat_eq(Astruct,[],c,[],param);
param.time = toc;
if ~isempty(param.S)
    [~, Rsx] = qr([Astruct.multiply(param.Xleft)*param.S,  param.Xleft, c], 'econ');
    [~, Rdx] = qr([param.Xleft, Astruct.multiply(param.Xleft)*(param.S)', c], 'econ');
    param.true_res = norm(Rsx * Rdx', 'fro')/norm(c)^2;
else
    [~, Rsx] = qr([Astruct.multiply(param.Xleft),  param.Xleft, c], 'econ');
    [~, Rdx] = qr([param.Xleft, Astruct.multiply(param.Xleft), c], 'econ');
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
L = ichol(-A);
ek_param.precond_left = -L;                                  % L\(A*x) = L\b
ek_param.precond_right = -L';                           % A*(U\y) = b, x = U\y
Astruct = struct_make(A, ek_param);
param = extended_ksm_mat_eq(Astruct,[],c,[],param);
param.time = toc;
if ~isempty(param.S)
    [~, Rsx] = qr([Astruct.multiply(param.Xleft)*param.S,  param.Xleft, c], 'econ');
    [~, Rdx] = qr([param.Xleft, Astruct.multiply(param.Xleft)*(param.S)', c], 'econ');
    param.true_res = norm(Rsx * Rdx', 'fro')/norm(c)^2;
else
    [~, Rsx] = qr([Astruct.multiply(param.Xleft),  param.Xleft, c], 'econ');
    [~, Rdx] = qr([param.Xleft, Astruct.multiply(param.Xleft), c], 'econ');
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
param = extended_ksm_mat_eq(A,[],c,[],param);
param.time = toc;
if ~isempty(param.S)
    [~, Rsx] = qr([A*(param.Xleft)*param.S,  param.Xleft, c], 'econ');
    [~, Rdx] = qr([param.Xleft, A*(param.Xleft)*(param.S)', c], 'econ');
    param.true_res = norm(Rsx * Rdx', 'fro')/norm(c)^2;
else
    [~, Rsx] = qr([A*(param.Xleft),  param.Xleft, c], 'econ');
    [~, Rdx] = qr([param.Xleft, A*(param.Xleft), c], 'econ');
    param.true_res = norm(Rsx * Rdx', 'fro')/norm(c)^2;
end
result_extended_backslash = param;
%%%
%%%%%




%%% restart
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
%%%


% append results
res.name = "4DLap";
res.size = j^2;
res.compress =  result_compress;
res.two_pass =  result_two_pass;
res.ratkrylov_backslash =  result_ratkrylov_backslash;
res.ratkrylov_plaincg =  result_ratkrylov_plaincg;
res.ratkrylov_icholcg =  result_ratkrylov_icholcg;
res.ratkrylov_multishiftcg =  result_ratkrylov_multishiftcg;
res.ratkrylov_multishiftcg_rem =  result_ratkrylov_multishiftcg_rem;
res.extended_plaincg =  result_extended_plaincg;
res.extended_icholcg =  result_extended_icholcg;
res.extended_backslash =  result_extended_backslash;
res.rest =  result_rest;
results{end+1} = res;
end
