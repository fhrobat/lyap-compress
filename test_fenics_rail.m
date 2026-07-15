%%
% Third experiment in the paper A. Casulli, F. Hrobat, D. Kressner "Lanczos with
% compression for symmetric Lyapunov equations": solving the generalized 
% symmetric Lyapunov equation A*X*E' + E*X*A' = -c*c' for FEniCS_Rail benchmark.
%%

T = table();

% we vary the max size of the mesh <-> vary size of A and E
for j = [79841, 317377, 1265537]

addpath(genpath('./'))

%generate A,E,L
% mult = @(v) -L\(A*((L')\v)) refers to the transformed Lyapunov

load('./matrices/rail_' + string(j)  + '.mat');

% A is negative definite
E = spdiags(sum(E, 2),0,j,j); 
L = chol(E, 'lower');
mult = @(v) -L\(A*((L')\v));
c = C(1,:)';

tol = 1e-4;
options.maxmem = 120;
options.tol = tol;
options.true_res = 0;



n = size(A,1);

% right-hand side after transformation
vec = L\c;

[Q, v1, v2, H] = full_orth_Arnoldi(mult, options.maxmem-1, vec/norm(vec));
alpha = diag(H);
beta = diag(H,-1);
options.eigmin = eigs(spdiags([[beta(1:end-1);0],alpha,[0;beta(1:end-1)]], [-1,0,1], length(alpha), length(alpha)), 1, 'smallestabs');
options.eigmax = eigs(spdiags([[beta(1:end-1);0],alpha,[0;beta(1:end-1)]], [-1,0,1], length(alpha), length(alpha)), 1, 'largestabs');
options.eigmin = options.eigmin /10;
options.eigmax = options.eigmax * 1.1;

%% compress

tic
[V, Y, result_compress] = lyap_compress(mult,vec,options);
result_compress.time = toc;

[~, R] = qr([mult(V) * Y, V, vec], 'econ');
n2 = size(Y, 2);
n1 = size(V, 2);
R1 = R(:, 1:n2);
R2 = R(:, n2+1:n2+n1);
R3 = R(:, end);
Rperm = [R2, R1, -R3];
result_compress.true_res = norm(R * Rperm', 'fro') / norm(vec)^2;


row = table();
row.N         = n;
row.method    = "compress";
row.time      = result_compress.time;
row.matvecs   = result_compress.iter;
row.true_res  = result_compress.true_res;
row.num_poles = result_compress.k;
T = [T; row];
clear V;
%%

%% two-pass
tic
[V, Y, result_two_pass] = two_pass_lanc(mult,vec,options);
result_two_pass.time = toc;

[~, R] = qr([mult(V) * Y, V, vec ], 'econ');
n2 = size(Y, 2);
n1 = size(V, 2);
R1 = R(:, 1:n2);
R2 = R(:, n2+1:n2+n1);
R3 = R(:, end);
Rperm = [R2, R1, -R3];
result_two_pass.true_res = norm(R * Rperm', 'fro') / norm(vec)^2;



row = table();
row.N         = n;
row.method    = "two-pass";
row.time      = result_two_pass.time;
row.matvecs   = result_two_pass.iter;
row.true_res  = result_two_pass.true_res;
row.num_poles = result_two_pass.k;
T = [T; row];
clear V;
%%




%% ratkrylov cycle

for num = [4,8,12]
xi = poles_ellipke(options.eigmin, options.eigmax, [], num);

for solver = ["plain", "ichol"]
options.precond = solver;
options.maxit = 10000;
options.tol_solver = 1e-8;

tic;
[V, Y, matvecs, n_poles] = lyap_ratkrylov_cycle(A, E, L, c, -xi(end:-1:1), options.tol, 30, options);
t_cycle = toc;

[~, R] = qr([mult(V) * Y, V, vec], 'econ');
n2 = size(Y, 2);
n1 = size(V, 2);
R1 = R(:, 1:n2);
R2 = R(:, n2+1:n2+n1);
R3 = R(:, end);
Rperm = [R2, R1, -R3];
res_cycle = norm(R * Rperm', 'fro') / norm(vec)^2;


row = table();
row.N         = n;
row.method    = "cycle_" + solver + "_" + string(num) +"_poles";
row.time      = t_cycle;
row.matvecs   =  matvecs;
row.true_res  = res_cycle;
row.num_poles = n_poles;
T = [T; row];
clear V;
end

end

%%

end