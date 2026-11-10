%%
% Sensitivity of compress to the approximation of the smallest eigenvalue for
% the 2D Laplacian of dimension N = 424^2 (tol = 1e-6). The smallest eigenvalue
% is multiplied by i = 1,...,10 and the Zolotarev poles are computed either with
% the number of poles recomputed from the perturbed interval (recomputed_k), or
% with the number of poles obtained from the exact interval (fixed_k).
%%

% set paths and folder for the results
root = setup_paths();
results_dir = get_results_dir(root);

% same problem as in test_2DLap.m
j = floor(10^2*sqrt(18));
A = (j+1)^2*gallery('tridiag',j);
A = kron(eye(j),A) + kron(A,eye(j));
c = zeros(j,1);
for i = 1:j
   c(i) = 1/(0.5*sqrt(2*pi))*exp(-(i/(j+1)-1/2)^2/(2*0.5^2));
end
c = kron(c,c);
eigmin = (j+1)^2 * (4 - 4 * cos(pi / (j+1)));
eigmax = (j+1)^2 * (4 - 4 * cos(j * pi/ (j+1)));

% normalize by norm(c)
normc = norm(c);
c = c/normc;
A = A / (normc^2);
mult = @(v) A*v;
eigmin = eigmin / (normc^2);
eigmax = eigmax / (normc^2);
fprintf("Sensitivity experiment for 2D Laplacian of dimension %d\n", j^2);

tol = 1e-6;
maxmem = 120;

% number of poles obtained with the exact spectral interval
k_exact = length(poles_Zolotarev(eigmin, eigmax, tol/(2*eigmax/eigmin)));

T = table();
for i = 1:10
    for variant = ["recomputed_k", "fixed_k"]
        options = [];
        options.tol = tol;
        options.maxmem = maxmem;
        options.eigmin = i*eigmin;
        options.eigmax = eigmax;
        if variant == "fixed_k"
            options.xi = poles_ellipke(options.eigmin, options.eigmax, [], k_exact);
        end

        [U, X, result] = lyap_compress(mult, c, options);

        % a posteriori scaled residual norm
        [~, R] = qr([mult(U) * X,  U, c], 'econ');
        n2 = size(X,2);
        n1 = size(U,2);
        R1 = R(:, 1:n2);
        R2 = R(:, n2+1:n2+n1);
        R3 = R(:, end);
        Rperm = [R2, R1, -R3];

        row = table();
        row.i         = i;
        row.variant   = variant;
        row.num_poles = result.k;
        row.matvecs   = result.iter;
        row.true_res  = norm(R * Rperm', 'fro')/norm(c)^2;
        T = [T; row];
        save_results(T, results_dir, 'test_2DLap_sensitivity');
        clear U;
    end
    disp(T(end-1:end, :));
end

disp(T);
