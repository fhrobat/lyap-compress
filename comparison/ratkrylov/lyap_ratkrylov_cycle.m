function [U, Y, matvecs, n_poles] = lyap_ratkrylov_cycle(A, E, L, c, poles, tol, m_max, options)
    % lyap_ratkrylov_cycle: Rational Krylov subspace method for the generalized
    % Lyapunov equation  A*X*E + E*X*A' = -c*c'.
    %
    % The equation is reduced to standard form via the substitution
    %   X_tilde = L'*X*L,   c_tilde = L\c,   A_tilde = L\A*L^{-T}
    % where E = LL' is the Cholesky decomposition
    %
    %
    % Inputs:
    %   A        : n x n Symmetric Negative Definite matrix
    %   E        : n x n Symmetric Positive Definite matrix
    %   L        : n x n lower Cholesky factor of E
    %   c        : n x 1 right-hand side vector
    %   poles    : array of real positive shift parameters (cycled)
    %   tol      : stopping tolerance on the relative residual norm
    %   m_max : maximum number of outer cycles
    %   options  : struct with fields
    %                .precond    – 'direct' | 'plain' | 'ichol'
    %                .tol_solver – inner iterative solver tolerance
    %                .maxit      – inner iterative solver max iterations
    %
    % Outputs:
    %   U          : n x m 
    %   Y          : m x m  solution of the projected Lyapunov equation
    %   matvecs    : total inner-solver iterations accumulated
    %   n_poles    : total number of poles used

    n = size(A, 1);
    num_poles = length(poles);

    c_tilde = L\c;
    norm_c = norm(c_tilde);

    max_steps = m_max * num_poles;
    U = zeros(n, max_steps + 1);
    H = zeros(max_steps + 1, max_steps);

    U(:, 1) = c_tilde/norm_c;

    step = 0;
    matvecs = 0;

    % Cache for factorizations / preconditioners
    R_cache = cell(num_poles, 1);
    Lichol_cache = cell(num_poles, 1);

    A_v_current = L\(A*(L'\U(:, 1)));

    % Projected matrix T_m
    Tm = [];

    for m = 1:m_max

        %% expand rational Krylov using each pole once
        for p = 1:num_poles
            step = step + 1;
            s_curr = poles(p);
            u_curr = U(:, step);
            switch options.precond

                case 'direct'
                    if m == 1
                        % Build and cache the Cholesky decomposition of (-A + s*E)
                        
                        R_cache{p} = decomposition(-A + s_curr * E, 'chol');
                    end
                    w_next = L' * (R_cache{p} \ (-L * u_curr));

                case 'plain'
                    % Unpreconditioned PCG 
                    M =  -A + s_curr * E;
                    b = -L * u_curr;
                    [w_next, ~, ~, it] = pcg(M, b, options.tol_solver, options.maxit);
                    w_next  = L' * w_next;
                    matvecs = matvecs + it;

                case 'ichol'
                    % if options.reorder == 1
                    %     if m == 1 
                    %         perm = dissect(A);
                    %         A_perm = A(perm,perm);
                    %         E_perm = E(perm,perm);
                    %     end
                    %     M_shifted = -A_perm + s_curr * E_perm;
                    %     if m == 1
                    %         icholopt.type   = 'ict';
                    %         icholopt.droptol = 1e-3;
                    %         Lichol_cache{p} = ichol(M_shifted, icholopt);
                    %     end
                    %     Lichol = Lichol_cache{p};
                    %     [w_next, ~, ~, it] = pcg(M_shifted, -L * u_curr(perm), ...
                    %                              options.tol_solver, options.maxit, ...
                    %                              Lichol, Lichol');
                    %     w_next(perm) = w_next;
                    %     w_next  = L' * w_next;
                    %     matvecs = matvecs + it;
                    % else
                        M = -A + s_curr * E;
                        if m == 1
                            icholopt.type   = 'ict';
                            icholopt.droptol = 1e-3;
                            Lichol_cache{p} = ichol(M ,icholopt);
                        end
                        Lichol = Lichol_cache{p};
                        b = -L * u_curr;
                        [w_next, ~, ~, it] = pcg(M, b, options.tol_solver, options.maxit, Lichol, Lichol');
                        w_next  = L' * w_next;
                        matvecs = matvecs + it;
                    % end
            end


            % CGS2
            h_col = U(:, 1:step)' * w_next;
            w_next = w_next - U(:, 1:step) * h_col;
            h_corr = U(:, 1:step)' * w_next;
            H(1:step, step) = h_col + h_corr;
            w_next = w_next - U(:, 1:step) * h_corr;
            H(step+1, step) = norm(w_next);
            U(:, step+1) = w_next / H(step+1, step);
            h_next = H(step+1, step);
            u_next = U(:, step+1);


            % update of the projected matrix T_m
            t_col = U(:, 1:step)' * A_v_current;

            if step == 1
                Tm = t_col;
            else
                Tm = [Tm,            t_col(1:end-1); ...
                      t_col(1:end-1)', t_col(end)];
            end
            Tm = (Tm + Tm') / 2;  


            A_v_next = L \ (A * (L' \ u_next));
            A_v_current = A_v_next;   
        end
        %%

        % solve projected equation
        Um = U(:, 1:step);

        e1 = zeros(step, 1);
        e1(1) = 1;
        bb_proj = norm_c^2 * (e1 * e1');

        Ym = lyap(Tm, bb_proj);

        %residual estimation
        Hm = H(1:step, 1:step);

        em = zeros(step, 1);
        em(step)  = 1;
        g = A_v_next - Um * (Um' * A_v_next);
        u2 = s_curr * u_next - g;
        u1_coeff = (Hm') \ (em * h_next);
        u1 = Um * (Ym * u1_coeff);

        [~, S_mat] = qr([u1, u2], 0);
        J = [0, 1; 1, 0];
        current_res = norm(S_mat * J * S_mat', 'fro') / norm_c^2;


        if current_res < tol
            break;
        end

    end  
    U = U(:, 1:step);
    Y = Ym;
    n_poles = step;
end