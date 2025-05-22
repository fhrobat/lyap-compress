function param = extended_ksm_lyap(A, C, param)
% param = EXTENDED_KSM_LYAP(A, C, param) approximates the solution of the
% Lyapunov equation
%
%                       A * X + X * A' + C * C' = 0
%
% via the extended Krylov method with the parameters specified in the param
% struct.
%
% A detailed description of the variables and parameter specifications can
% be found in EXTENDED_KSM_MAT_EQ.
%
% See also EXTENDED_KSM_SYLV and EXTENDED_KRYLOV.

%%
% Initialize
Acount = 0;
flag = 0;
matvecs = 0;
residual_norms = [];

res_init = norm(C'*C,param.norm);                                           % scale by initial residual

% Solution storage
r = size(C, 2);
Xleft = [];
S  = [];
Cold = C;

% Set outer memory allowance. Note that EXTENDED_KRYLOV expands the space 2
% block vectors at a time.
max_iterations = floor(param.memory_max/2/r);
if max_iterations == 0
    error('Ran out of memory. Increase and try again.')
end

% Initialize inner solver memory allowance
param.max_iterations = floor(param.memory_max/r);
switch A.solver
    case {'bgmres','gmres'}
        % Have to limit number of basis vectors we can generate before
        % restarting.
        if param.max_iterations <= 0
            error('Ran out of memory for A solves. Increase and try again.')
        end

    case {'bcg', 'cg'}
        % Just have to make sure we have enough storage for 3 vectors
        % at a time
        if param.max_iterations >= 3
            param.max_iterations = inf;
        else
            error('Ran out of memory for A solves. Increase and try again.')
        end
end

for j = 1:max_iterations
    % Augment the Krylov subspace, two block vectors at a time. Note that
    % EXTENDED_KRYLOV is written to reset Acount and matvecs each time it
    % is called, so we need not worry about redundancy.
    param = extended_krylov(A, C, param);
    Acount = Acount + param.Acount;
    matvecs = matvecs + param.matvecs;
    
    % Retrieve the projection of A in the Krylov subspaces
    Am = param.H/param.K(1:end - r, :);

    % Solve the projected matrix equation    
    Y = lyap(Am(1:end-r,:), ...
        (param.E1 * param.H0) * (param.E1 * param.H0)');

    % Safeguard to preserve symmetry
    Y = (Y + Y')/2;

    % Compute approximate residual
    switch param.norm
        case 'fro'
            space_dim = size(param.V,2);
            I = eye(space_dim - r);
            factor = [param.V(:, 1:space_dim - r)*...
                (Y*(param.K(1:space_dim - r, :)'\...
                I(:,end - r+1:end))),...
                -struct_mult(A,param.V(:, end - r+1:end))*...
                param.K(end - r+1:end, end - r+1:end)+...
                param.V(:, end - r+1:end)*...
                param.H(end - r+1:end, end - r+1:end)];
            [~,R] = qr(factor,0);
            J = zeros(2*r);
            J(1:r,r+1:2*r) = eye(r);
            J(r+1:2*r,1:r) = eye(r);
            res_approx = norm(R*J*R','fro')/res_init;
            
        case {2, '2'}
            space_dim = size(param.V,2);
            I = eye(space_dim - r);
            factor = [param.V(:, 1:space_dim - r)*...
                (Y*(param.K(1:space_dim - r, 1:space_dim - r)'\...
                I(:,end - r + 1:end))), ...
                -struct_mult(A, param.V(:, end - r+1:end))*...
                param.K(end - r+1:end, end - r+1:end)+...
                param.V(:, end - r+1:end)*...
                param.H(end - r+1:end, end - r+1:end)];
            [~,R] = qr(factor,0);
            J = zeros(2*r);
            J(1:r,r+1:2*r) = eye(r);
            J(r+1:2*r,1:r) = eye(r);
            res_approx = norm(R*J*R')/res_init; 
    end

    % Decompose Y for solution update
    [W, Snew] = eig(Y);

    % Print residual
    switch param.verbose
        case 0   % Don't print anything
            residual_norms = [residual_norms, res_approx];

        case 1   % Print the approximate residual
            fprintf('     Inner iter: %d, res = %e\n', j, res_approx);
            residual_norms = [residual_norms, res_approx];

        case 2   % Calculate and print the true residual
            YY = W * diag(sqrt(diag(Snew)));

            if j == 1
                Xlefttemp = param.V(:, 1:end - r) * YY;
            else
                Xlefttemp = [Xleft, param.V(:, 1:end - r) * YY];
            end

            if ~isempty(S)
                fact = [struct_mult(A,Xlefttemp)*...
                    blkdiag(S, eye(size(Xlefttemp, 2) -size(S, 1))),...
                    Xlefttemp, Cold];
            else
                fact = [struct_mult(A,Xlefttemp), Xlefttemp, Cold];
            end

            Theta = eye(size(fact, 2));
            Theta(1:2 * size(Xlefttemp, 2), 1:2 * size(Xlefttemp,2))= ...
                [zeros(size(Xlefttemp, 2)), eye(size(Xlefttemp, 2));...
                eye(size(Xlefttemp, 2)), zeros(size(Xlefttemp, 2))];
            [~, R] = qr(full(fact), 0);

            switch param.norm
                case 'fro'
                    res_true = norm(full(R * Theta * R'), 'fro')/res_init;

                case {2, '2'}
                    res_true = norm(full(R * Theta * R'))/res_init;
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
    
    % Reduce allowable basis size for the inner solver by 2, since vectors
    % corresponding to both inf and 0 poles are computed per basis-building
    % call
    param.max_iterations = param.max_iterations - 2;
    switch A.solver
        case {'bgmres','gmres'}
            if param.max_iterations <= 0
                error('Ran out of memory for A solves. Increase and try again.')
            end

        case {'bcg', 'cg'}
            if param.max_iterations >= 3
                param.max_iterations = inf;
            else
                error('Ran out of memory for A solves. Increase and try again.')
            end
    end
end
if res_approx >= param.tol_res && j == max_iterations
    flag = 2;
end

% Update solution
Xleft = [Xleft, param.V(:,1:end-r) * W];

% Compress solution
[Xleft, S] = symm_compress_factors(Xleft, blkdiag(S, Snew),...
    param);

if min(eig(S)) < 0   % Clean the solution if it is indefinite
    [Xleft, S] = symm_compress_clean(Xleft, S);
end


% Update param outputs alphabetically
param.Acount = Acount;
param.flag = flag;
param.matvecs = matvecs;
param.residual_norms = residual_norms;
param.S = S;
param.Xleft = Xleft;
end
