function [V, Y, shifts, matvecs] = lyap_ratkrylov_adaptive(A, E, L, c, m_max, options)
    % lyap_ratkrylov_adaptive: Rational Krylov subspace method for the generalized
    % Lyapunov equation  A*X*E + E*X*A' = -c*c'
    % with ADAPTIVE pole selection.
    %
    % The equation is reduced to standard form via the substitution
    %   X_tilde = L'*X*L,   c_tilde = L\c,   A_tilde = L\A*L^{-T}
    % where E = LL' is the Cholesky decomposition
    %
    %
    % Inputs:
    %   A       : n x n Symmetric Negative Definite matrix
    %   E       : n x n Symmetric Positive Definite matrix
    %   L       : n x n lower Cholesky factor of E
    %   c       : n x 1 right-hand side vector
    %   m_max   : maximum allowed Krylov subspace dimension
    %   options : struct with fields
    %               .eigmin     – lower bound for the mirrored spectrum  (positive scalar)
    %               .eigmax     – upper bound for the mirrored spectrum  (positive scalar)
    %               .tol        – stopping tolerance on the relative residual norm
    %               .precond    – 'direct' | 'plain' | 'ichol'
    %               .tol_solver – inner iterative solver tolerance
    %               .maxit      – inner iterative solver max iterations
    %
    % Outputs:
    %   V      : n x m 
    %   Y      : m x m  solution of the projected Lyapunov equation
    %   shifts : (m+1) x 1  poles used
    %   matvecs: total inner-solver iterations accumulated

    n = size(A, 1);
    s0_1 = options.eigmin;   % initial lower bound of the mirrored spectrum
    s0_2 = options.eigmax;   % initial upper bound of the mirrored spectrum
    tol = options.tol;
    matvecs = 0;


    V = zeros(n, m_max + 1);
    H = zeros(m_max + 1, m_max);


    c_tilde = L\c;
    norm_c  = norm(c_tilde);
    V(:, 1) = c_tilde / norm_c; %first pole is chosen as infinity

    shifts = [];


    A_v_current = L\(A*(L'\V(:, 1)));
    t_col = V(:,1)'*A_v_current;

    % Projected matrix T_m = Vm' * A_tilde * Vm
    T_m = [];


    for m = 1:m_max
        Vm = V(:,1:m);

        % Append one column/row to T_{m-1} to obtain T_m.

        if m == 1
            T_m = t_col;
        else
            T_m = [T_m,             t_col(1:m-1); ...
                   t_col(1:m-1)',   t_col(m)];
        end
        T_m = (T_m+T_m')/2; 


        % Solve projected Lyapunov equation
        c_m = zeros(m, 1);
        c_m(1) = norm_c;
        rhs_m = -(c_m * c_m');
        Y = lyap(T_m, rhs_m);

        % Use the Ritz values of T_m to update the spectral bounds and compute
        % the next pole.
        if m == 1
            s_next = s0_1;
        else
            lambdas = eig(T_m, 'vector');        % real, strictly negative
            s0_1 = min(s0_1, min(-lambdas));    % tighten lower spectral bound
            s0_2 = max(s0_2, max(-lambdas));    % tighten upper spectral bound
            s_next  = select_adaptive_symm(lambdas, shifts, s0_1, s0_2);
        end
        shifts = [shifts; s_next];

        % extend rational Krylov subspace
        if strcmp(options.precond, 'direct')
            % Direct sparse solve
            w = L'*((A - s_next * E)\(L * V(:, m)));

        elseif strcmp(options.precond, 'plain')
            % Unpreconditioned PCG
            [w, ~, ~, it] = pcg(-A + s_next * E, -L * V(:, m), options.tol_solver, options.maxit);
            w = L' * w;
            matvecs = matvecs + it;

        elseif strcmp(options.precond, 'ichol')
            % PCG with incomplete Cholesky preconditioner
            M = -A + s_next * E;
            b = -L*V(:, m);
            icholopt.type = 'ict';
            icholopt.droptol = 1e-3;
            Lichol = ichol(M, icholopt);
            [w, ~, ~, it] = pcg(M, b, options.tol_solver, options.maxit, Lichol, Lichol');
            w       = L' * w;
            matvecs = matvecs + it;
        end

        % CGS2
        h_tmp = V(:, 1:m)' * w;
        H(1:m, m) = h_tmp;
        w = w - V(:, 1:m) * h_tmp;
        h_corr = V(:, 1:m)' * w;
        H(1:m, m) = H(1:m, m) + h_corr;
        w = w - V(:, 1:m) * h_corr;

        %extend H and V
        H(m+1, m) = norm(w);
        V(:, m+1) = w / H(m+1, m);

        A_v_next = L \ (A * (L' \ V(:, m+1)));
        
        v_proj = Vm' * A_v_next;
       

        % next column of T_m
        t_col = [v_proj; V(:,m+1)' * A_v_next];

        % redisual estimation
        e_m_col = zeros(m, 1);
        e_m_col(m) = 1;

        
        z_m = Y * (H(1:m, 1:m)' \ e_m_col);
        u1  = Vm * z_m * H(m+1, m);
        g = A_v_next - Vm * v_proj;
        u2 = V(:, m+1) * s_next - g;
        [~, S_qr] = qr([u1, u2], 0);
        J = [0, 1; 1, 0]; 
        res_norm = norm(S_qr * J * S_qr', 'fro');

        if res_norm / norm_c^2 < tol
            break;
        end

    end

    V = V(:, 1:m);
    Y = -Y;
end



% AUXILIARY FUNCTION  –  adaptive pole selection
function s_new = select_adaptive_symm(lambdas, current_shifts, s0_1, s0_2)

    test_points = logspace(log10(s0_1), log10(s0_2), 500);
    num = prod(abs(test_points - current_shifts), 1);
    den = prod(abs(test_points - (lambdas)),     1);
    vals = num ./ den;

    [~, idx] = max(vals);
    s_new = test_points(idx);
end