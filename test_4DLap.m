%%
% First experiment in the paper A. Casulli, F. Hrobat, D. Kressner "Lanczos with
% compression for symmetric Lyapunov equations": solving the symmetric 
% Lyapunov equation A*X + X*A' = -c*c', where A is the 2D Laplacian and 
% c is a 2D Gaussian vector.
%%

% table that contains all the results
T = table();

% j is the square root of N

for j= [floor(10^2*sqrt(18)),floor(10^2 * sqrt(36)), floor(10^2 * sqrt(72)), 1200] 

addpath(genpath('./'))

% generate function handle for mult, right-hand side and eigenvalues
%%
A = (j+1)^2*gallery('tridiag',j);
A = kron(eye(j),A) + kron(A,eye(j));
mult = @(v) A*v;
c = zeros(j,1);
for i = 1:j
   c(i) = 1/(0.5*sqrt(2*pi))*exp(-(i/(j+1)-1/2)^2/(2*0.5^2));
end
c = kron(c,c);
options.eigmin = (j+1)^2 * (4 - 4 * cos(pi / (j+1)));
options.eigmax = (j+1)^2 * (4 - 4 * cos(j * pi/ (j+1)));
%%

% normalize by norm(c)
normc = norm(c);
c = c/normc;
mult = @(v) mult(v) / (normc)^2;
options.eigmin = options.eigmin / (normc^2);
options.eigmax = options.eigmax / (normc^2);
A = A / (normc^2);
fprintf("Experiments for 4D Laplacian of dimension %d\n", j^2);

% set options
tol = 1e-6;
options.maxmem = 120;
options.tol = tol;
options.true_res = 0;

%% compress
tic
[U, X, result_compress] = lyap_compress(mult,c,options);
result_compress.time = toc;

% true residual norm for compress
[~, R] = qr([mult(U) * X,  U, c], 'econ');
n2 = size(X,2);
n1 = size(U,2);
R1 = R(:, 1:n2);
R2 = R(:, n2+1:n2+n1);
R3 = R(:, end);
Rperm = [R2, R1, -R3];
result_compress.true_res = norm(R * Rperm', 'fro')/norm(c)^2;

row = table();
row.N         = j^2;
row.method    = "compress";
row.time      = result_compress.time;
row.matvecs   = result_compress.iter;
row.true_res  = result_compress.true_res;
row.num_poles = result_compress.k;
T = [T; row];
clear U;
%%


%% two-pass
tic
[U, X, result_two_pass] = two_pass_lanc(mult,c,options);
result_two_pass.time = toc;

% true residual norm for two-pass
[~, R] = qr([mult(U) * X,  U, c], 'econ');
[n1, n2] = size(X);
R1 = R(:, 1:n2);
R2 = R(:, n2+1:n2+n1);
R3 = R(:, end);
Rperm = [R2, R1, -R3];
result_two_pass.true_res = norm(R * Rperm', 'fro')/norm(c)^2;

row = table();
row.N         = j^2;
row.method    = "two-pass";
row.time      = result_two_pass.time;
row.matvecs   =  result_two_pass.iter;
row.true_res  = result_two_pass.true_res;
row.num_poles = result_two_pass.k;
T = [T; row];
clear U;
%%



%% restart
M.multiply = @(v) -mult(v);

param = [];
param.max_restarts = 11000/ options.maxmem;
param.memory_max = options.maxmem; % set same memory costraint
param.norm = 'fro';
param.tol_res = tol;
param.tol_comp = tol*1e-4; % default tolerance
param.verbose = 0;

tic
param = restarted_ksm_mat_eq(M,[],c,[],param);
param.time_rest = toc;

% true residual norm for restart
if ~isempty(param.S)
    [~, Rsx] = qr([M.multiply(param.Xleft)*param.S,  param.Xleft, c], 'econ');
    [~, Rdx] = qr([param.Xleft, M.multiply(param.Xleft)*(param.S)', c], 'econ');
    param.true_res = norm(Rsx * Rdx', 'fro')/norm(c)^2;
else
    [~, Rsx] = qr([M.multiply(param.Xleft),  param.Xleft, c], 'econ');
    [~, Rdx] = qr([param.Xleft, M.multiply(param.Xleft), c], 'econ');
    param.true_res = norm(Rsx * Rdx', 'fro')/norm(c)^2;
end


% result for restart

row = table();
row.N         = j^2;
row.method    = "restart";
row.time      = param.time_rest;
row.matvecs   =  param.matvecs;
row.true_res  = param.true_res;
row.num_poles = NaN;
T = [T; row];
clear param;

%%


%% ratkrylov cycle
xi = poles_ellipke(options.eigmin, options.eigmax, [], 4);

for solver = ["direct", "plain", "ichol"]
options.precond = solver;
options.maxit = 10000;
options.tol_solver = 1e-10;

tic;
[V, Y, matvecs, n_poles] = lyap_ratkrylov_cycle(-A, speye(j^2), speye(j^2), c, -xi(end:-1:1), options.tol, 30, options);
t_cycle = toc;

[~, R] = qr([mult(V) * Y, V, c], 'econ');
n2 = size(Y, 2);
n1 = size(V, 2);
R1 = R(:, 1:n2);
R2 = R(:, n2+1:n2+n1);
R3 = R(:, end);
Rperm = [R2, R1, -R3];
res_cycle = norm(R * Rperm', 'fro') / norm(c)^2;


row = table();
row.N         = j^2;
row.method    = "cycle_" + solver;
row.time      = t_cycle;
row.matvecs   =  matvecs;
row.true_res  = res_cycle;
row.num_poles = n_poles;
T = [T; row];
clear V;
end

%%

%% ratkrylov adaptive 

for solver = ["direct", "plain", "ichol"]

options.precond = solver;
options.maxit = 10000;
options.tol_solver = 1e-10;

tic;
[V, Y, shifts, matvecs] = lyap_ratkrylov_adaptive(-A, speye(j^2), speye(j^2), c, 70, options);
t_ada = toc;

[~, R] = qr([mult(V) * Y, V, c], 'econ');
n2 = size(Y, 2);
n1 = size(V, 2);
R1 = R(:, 1:n2);
R2 = R(:, n2+1:n2+n1);
R3 = R(:, end);
Rperm = [R2, R1, -R3];
res_ada = norm(R * Rperm', 'fro') / norm(c)^2;


row = table();
row.N         = j^2;
row.method    = "adaptive_" + solver;
row.time      = t_ada;
row.matvecs   =  matvecs;
row.true_res  = res_ada;
row.num_poles = n_poles;
T = [T; row];
clear V;
end

%% 

%% ratkrylov mutlishiftcg

options.maxit = 10000;
options.tol_solver = 1e-8;
options.variant = 'multishiftcg';

tic
c0 = options.eigmax/options.eigmin;
options.xi = poles_Zolotarev(options.eigmin, options.eigmax, options.tol/(c0 * 2));
[U, X, nmatvecs, k] =  lyap_ratkrylov_multishift(A,c,options);
result_ratkrylov_multishiftcg.time = toc;

% true residual norm for ratkrylov mutltishiftcg
[~, R] = qr([mult(U) * X,  U, c], 'econ');
[n1, n2] = size(X);
R1 = R(:, 1:n2);
R2 = R(:, n2+1:n2+n1);
R3 = R(:, end);
Rperm = [R2, R1, -R3];
result_ratkrylov_multishiftcg.true_res = norm(R * Rperm', 'fro')/norm(c)^2;


row = table();
row.N         = j^2;
row.method    = "multishiftCG";
row.time      = result_ratkrylov_multishiftcg.time;
row.matvecs   =  nmatvecs;
row.true_res  = result_ratkrylov_multishiftcg.true_res;
row.num_poles = k;
T = [T; row];
clear U;
%%


%% ratkrylov mutlishiftcg with removal
options.variant = 'multishiftcg_rem';
tic
c0 = options.eigmax/options.eigmin;
options.xi = poles_Zolotarev(options.eigmin, options.eigmax, options.tol/(c0 * 2));
[U, X, nmatvecs, k] =  lyap_ratkrylov_multishift(A,c,options);
result_ratkrylov_multishiftcg_rem.time = toc;

% true residual norm for ratkrylov mutlishiftcg with removal
[~, R] = qr([mult(U) * X,  U, c], 'econ');
[n1, n2] = size(X);
R1 = R(:, 1:n2);
R2 = R(:, n2+1:n2+n1);
R3 = R(:, end);
Rperm = [R2, R1, -R3];
result_ratkrylov_multishiftcg_rem.true_res = norm(R * Rperm', 'fro')/norm(c)^2;



row = table();
row.N         = j^2;
row.method    = "multishiftCG_rem";
row.time      = result_ratkrylov_multishiftcg_rem.time;
row.matvecs   =  nmatvecs;
row.true_res  = result_ratkrylov_multishiftcg_rem.true_res;
row.num_poles = k;
T = [T; row];
clear U;
%%


%% extended
ek_param = [];
ek_param.solver = 'cg';
ek_param.symmetric = true;
ek_param.max_iterations = 10000;
ek_param.max_restarts = param.max_restarts;
ek_param.verbose = 0;
ek_param.tol_res = tol*1e-4;

param = [];
param.max_restarts = 11000/ options.maxmem;
param.memory_max = options.maxmem; % set same memory costraint
param.norm = 'fro';
param.tol_res = tol;
param.tol_comp = tol*1e-4; % default tolerance
param.verbose = 0;


%%% icholPCG
tic
icholopt.type   = 'ict';
icholopt.droptol = 1e-3;
L = ichol(A,icholopt);
ek_param.precond_left = -L;
ek_param.precond_right = -L';
Astruct = struct_make(-A, ek_param);
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

row = table();
row.N         = j^2;
row.method    = "extended-icholPCG";
row.time      = param.time;
row.matvecs   =  param.matvecs;
row.true_res  = param.true_res;
row.num_poles = param.Acount/2;
T = [T; row];
%%%


param = [];
param.max_restarts = 11000/ options.maxmem;
param.memory_max = options.maxmem; % set same memory costraint
param.norm = 'fro';
param.tol_res = tol;
param.tol_comp = tol*1e-4; % default tolerance
param.verbose = 0;


%%% backslash
tic
param = extended_ksm_mat_eq(-A,[],c,[],param);
param.time = toc;
if ~isempty(param.S)
    [~, Rsx] = qr([-A*(param.Xleft)*param.S,  param.Xleft, c], 'econ');
    [~, Rdx] = qr([param.Xleft, -A*(param.Xleft)*(param.S)', c], 'econ');
    param.true_res = norm(Rsx * Rdx', 'fro')/norm(c)^2;
else
    [~, Rsx] = qr([-A*(param.Xleft),  param.Xleft, c], 'econ');
    [~, Rdx] = qr([param.Xleft, -A*(param.Xleft), c], 'econ');
    param.true_res = norm(Rsx * Rdx', 'fro')/norm(c)^2;
end
row = table();
row.N         = j^2;
row.method    = "extended-backslash";
row.time      = param.time;
row.matvecs   =  param.matvecs;
row.true_res  = param.true_res;
row.num_poles = param.Acount/2;
T = [T; row];
%%%
%%

end
