%%
% Second experiment in the paper A. Casulli, F. Hrobat, D. Kressner "Lanczos with
% compression for symmetric Lyapunov equations": solving the symmetric 
% Lyapunov equation A*X + X*A' = -c*c', where A is the 3D Laplacian and 
% c is a 3D Gaussian vector.
%%

% table that contains all the results
T = table();

% j is the 3rd root of N

for j= [floor(10^(4/3)*18^(1/3)),floor(10^(4/3)*36^(1/3)), floor(10^(4/3)*72^(1/3)), floor(10^(4/3)*144^(1/3))] 

addpath(genpath('./'))

% generate function handle for mult, right-hand side and eigenvalues

A = (j+1)^2 * gallery('tridiag', j); 


A = kron(kron(speye(j), speye(j)), A) + ...
    kron(kron(speye(j), A), speye(j)) + ...
    kron(kron(A, speye(j)), speye(j));

mult = @(v) A*v;

c = zeros(j,1);
for i = 1:j
    c(i) = 1/(0.5*sqrt(2*pi)) * exp(-(i/(j+1) - 1/2)^2 / (2*0.5^2));
end
c = kron(kron(c, c), c); 

options.eigmin = (j+1)^2 * (6 - 6*cos(pi/(j+1)));
options.eigmax = (j+1)^2 * (6 - 6*cos(j*pi/(j+1)));

% normalize by norm(c)
normc = norm(c);
c = c/normc;
mult = @(v) mult(v) / (normc)^2;
options.eigmin = options.eigmin / (normc^2);
options.eigmax = options.eigmax / (normc^2);
A = A / (normc^2);
fprintf("Experiments for 6D Laplacian of dimension %d\n", j^3);

% set options
tol = 1e-12;
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
row.N         = j^3;
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
row.N         = j^3;
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
param.tol_comp = eps;
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
row.N         = j^3;
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

for solver = [ "plain", "ichol"]
options.precond = solver;
options.maxit = 10000;
options.tol_solver = 1e-13;

tic;
[V, Y, matvecs, n_poles] = lyap_ratkrylov_cycle(-A, speye(j^3), speye(j^3), c, -xi(end:-1:1), options.tol, 30, options);
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
row.N         = j^3;
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

for solver = [ "plain", "ichol"]

options.precond = solver;
options.maxit = 10000;
options.tol_solver = 1e-13;

tic;
[V, Y, shifts, matvecs] = lyap_ratkrylov_adaptive(-A, speye(j^3), speye(j^3), c, 70, options);
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
row.N         = j^3;
row.method    = "adaptive_" + solver;
row.time      = t_ada;
row.matvecs   =  matvecs;
row.true_res  = res_ada;
row.num_poles = n_poles;
T = [T; row];
clear V;

end

%% 


%% extended
ek_param = [];
ek_param.solver = 'cg';
ek_param.symmetric = true;
ek_param.max_iterations = 10000;
ek_param.max_restarts = param.max_restarts;
ek_param.verbose = 0;
ek_param.tol_res = 1e-13;

param = [];
param.max_restarts = 11000/ options.maxmem;
param.memory_max = options.maxmem; % set same memory costraint
param.norm = 'fro';
param.tol_res = tol;
param.tol_comp = eps; 
param.verbose = 0;

%%% plainCG
tic
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
row.N         = j^3;
row.method    = "extended-plainCG";
row.time      = param.time;
row.matvecs   =  param.matvecs;
row.true_res  = param.true_res;
row.num_poles = param.Acount/2;
T = [T; row];
%%%




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
row.N         = j^3;
row.method    = "extended-icholPCG";
row.time      = param.time;
row.matvecs   =  param.matvecs;
row.true_res  = param.true_res;
row.num_poles = param.Acount/2;
T = [T; row];
clear param;
%%%

%%
end
