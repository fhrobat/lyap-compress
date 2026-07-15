function param = extended_ksm_sylv(A, B, C, D, param)
% param = EXTENDED_KSM_SYLV(A, B, C, D, param) approximates the solution
% of the Sylvester equation
%
%                       A * X + X * B + C * D' = 0
%
% via the extended Krylov method with the parameters specified in the param
% struct.
%
% A detailed description of the variables and parameter specifications can
% be found in EXTENDED_KSM_MAT_EQ.
% 
% See also EXTENDED_KSM_LYAP and EXTENDED_KRYLOV.

%%
% Initialize
Acount = [0; 0];
flag = 0;
matvecs = [0; 0];
residual_norms = [];

res_init = norm(C'*D,param.norm);

r = size(C, 2);                                                             % assume same as D
Xleft = []; 
Xright = [];

% Set up distinct param structs for each space.
paramA = param;
paramB = param;

% Set outer memory allowance; assume each space has the same memory
% limitations. Note that since EXTENDED_KRYLOV expands the space 2 block
% vectors at a time, we must consider that 4 basis vectors are generated
% each iteration, 2 for each space.
max_iterations = floor(param.memory_max/4/r);
if max_iterations == 0
    error('Ran out of memory. Increase and try again.')
end

% Initialize inner solver memory allowances separately for A and B; we
% assume the spaces are generated sequentially, so that each has access to
% the same amount of memory per outer iteration
paramA.max_iterations = floor(param.memory_max/2/r);
switch A.solver
    case {'bgmres','gmres'}
        % Have to limit number of basis vectors we can generate before
        % restarting. Note that backslash doesn't use this parameter.
        if paramA.max_iterations <= 0
            error('Ran out of memory for A solves. Increase and try again.')
        end
    case {'bcg', 'cg'}
        % Just have to make sure we have enough storage for 3 vectors
        % at a time
        if paramA.max_iterations >= 3
            paramA.max_iterations = inf;
        else
            error('Ran out of memory for A solves. Increase and try again.')
        end
end
if strcmp(A.solver, 'backslash')
    paramB.max_iterations = paramA.max_iterations;
else
    paramB.max_iterations = floor(param.memory_max/2/r);
end
switch B.trans.solver
    case {'bgmres','gmres'}
        if paramB.max_iterations <= 0
            error('Ran out of memory for B solves. Increase and try again.')
        end
    case {'bcg', 'cg'}
        if paramB.max_iterations >= 3
            paramB.max_iterations = inf;
        else
            error('Ran out of memory for B solves. Increase and try again.')
        end
end

for j = 1:max_iterations
    % Augment the Krylov spaces
    paramA = extended_krylov(A, C, paramA);
    if isnumeric(B)
        paramB = extended_krylov(B', D, paramB);
    else
        paramB = extended_krylov(B.trans, D, paramB);
    end
    Acount = Acount + [paramA.Acount; paramB.Acount];
    matvecs = matvecs + [paramA.matvecs; paramB.matvecs];
    
    % Retrieve the projection of A in the Krylov subspaces
    Am = paramA.H/paramA.K(1:end - r, :);
    Bm = paramB.H/paramB.K(1:end - r, :);

    % Solve the projected matrix equation
    Y = lyap(Am(1:end-r,:), Bm(1:end-r,:)',...
            (paramA.E1 * paramA.H0) * (paramB.E1 * paramB.H0)');

    % Compute approximate residual
    switch param.norm
        case 'fro'
            res_approx = sqrt(...
                norm(Am(end - r + 1:end, :) * Y,'fro')^2 +...
                norm(Y * Bm(end - r + 1 : end, :)','fro')^2)...
                /res_init;

        case {2, '2'}
            res_approx = max(norm(Am(end - r + 1:end, :) * Y), ...
            norm(Y * Bm(end - r + 1 : end, :)'))/res_init;
    end

    % Print residual
    switch param.verbose
        case 0    % Don't print anything
            residual_norms = [residual_norms, res_approx];

        case 1  % Print the approximate residual
            fprintf('     Inner iter: %d, res = %e\n', j, res_approx);
            residual_norms = [residual_norms, res_approx];

        case 2   % Calculate and print the true residual
            X = paramA.V(:, 1:end-r) * Y * paramB.V(:, 1:end-r)';

            if isnumeric(B)
                res_true = norm(struct_mult(A,X) +...
                    (B' * X')' + C * D', param.norm)/...
                    res_init;

            elseif isstruct(B)
                res_true = norm(struct_mult(A,X) +...
                    struct_mult(B.trans, X')' + C * D', param.norm)/...
                    res_init;
            end
            fprintf('     Inner iter: %d, res = %e, true_res = %e\n',...
                j, res_approx, res_true);
            residual_norms = [residual_norms, res_true];
    end

    % Convergence check with approximate residual
    if res_approx < param.tol_res
        flag = 1;
        break
    end
    
    % Reduce allowable basis size for the inner solvers by 2, since vectors
    % corresponding to both inf and 0 poles are computed per basis-building
    % call
    paramA.max_iterations = paramA.max_iterations - 2;
    switch A.solver
        case {'bgmres','gmres'}
            if paramA.max_iterations <= 0
                error('Ran out of memory for A solves. Increase and try again.')
            end
        case {'bcg', 'cg'}
            if paramA.max_iterations >= 3
                paramA.max_iterations = inf;
            else
                error('Ran out of memory for A solves. Increase and try again.')
            end
    end
    if strcmp(A.solver, 'backslash')
        paramB.max_iterations = paramA.max_iterations;
    else
        paramB.max_iterations = paramB.max_iterations - 2;
    end
    switch B.trans.solver
        case {'bgmres','gmres'}
            if paramB.max_iterations <= 0
                error('Ran out of memory for B solves. Increase and try again.')
            end
        case {'bcg', 'cg'}
            if paramB.max_iterations >= 3
                paramB.max_iterations = inf;
            else
                error('Ran out of memory for B solves. Increase and try again.')
            end
    end
end
if res_approx >= param.tol_res && j == max_iterations
    flag = 2;
end

% Update the solution
[W, S, Z] = svd(Y);
Xleft = [Xleft, paramA.V(:,1:end-r) * W * sqrt(S)];
Xright = [Xright, paramB.V(:,1:end-r) * Z * sqrt(S)];

% Compress the solution
[Xleft, Xright] = compress_factors(Xleft, Xright, param);

% Update param outputs alphabetically
param.Acount = Acount;
param.flag = flag;
param.matvecs = matvecs;
param.residual_norms = residual_norms;
param.Xleft = Xleft;
param.Xright = Xright;
end