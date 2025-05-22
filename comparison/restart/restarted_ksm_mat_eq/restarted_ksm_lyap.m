function param = restarted_ksm_lyap(A, C, param)
% param = RESTARTED_KSM_LYAP(A, C, param) approximates the solution of the
% Lyapunov equation
%
%                       A * X + X * A' + C * C' = 0
%
% via the restarted polynomial Krylov method with the parameters specified
% in the param struct.
%
% A detailed description of the variables and parameter specifications can
% be found in RESTARTED_KSM_MAT_EQ.
%
% See also RESTARTED_KSM_SYLV and BLOCK_ARNOLDI.

%%
% Initialize
Acount = 0;
flag = 0;
matvecs = 0;
residual_ranks = [];
solution_ranks = [];
residual_norms = [];
it_per_cycle = [];

switch param.norm
        case 'fro'
            res_init = norm(C'*C,param.norm);                                           % scale by initial residual   
        case {2, '2'}
            [~,RC]=qr(C,0);
            res_init=norm(RC*RC');
end    
% Solution storage
r = size(C, 2);
Xleft = [];
S  = [];
L = eye(r);
Cold = C;
min_core_eig = [];

for k = 1:param.max_restarts
    if param.verbose == 1 || param.verbose == 2
        fprintf('Restart cycle: %d----------------------\n',k);
    end
    
    if k > 1    % Clear previous Krylov basis
        param = clear_basis(param);
    end
        
    if r == 0
        error('Compression tolerance is too high. Lower and try again.')
    end
    
    % Adjust memory allowance
    param.max_iterations = floor(param.memory_max/r);                       % Need to store in param in order to feed to block_arnoldi
    if param.max_iterations == 0
        error('Ran out of memory. Increase and try again.')
    end
    
    for j = 1:param.max_iterations
		% Augment the Krylov subspace
        param = block_arnoldi(A, C, param);
        Acount = Acount + 1;
        matvecs = matvecs + r;
        
        % Solve the projected matrix equation
        Y = lyap(param.Hshort,...
            param.E1 * param.H0 * L * param.H0' * param.E1');
        
      	% Safeguard to preserve symmetry
		Y = (Y + Y')/2;
        
        % Compute approximate residual
        switch param.norm
            case 'fro'
                res_approx = sqrt(2)*...
                    norm(param.Hnext*Y(end-r+1:end,:),'fro')/res_init;
            case {2, '2'}
                res_approx = norm(param.Hnext*Y(end-r+1:end,:))/res_init;
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
                YY = W ;

                if k == 1            		
                    Xlefttemp = param.Vshort * YY;
                else
                    Xlefttemp = [Xleft, param.Vshort * YY];
                end

                if ~isempty(S)
                    fact = [struct_mult(A,Xlefttemp)*...
                        blkdiag(S, Snew),...
                        Xlefttemp, Cold];
                else
                    fact = [struct_mult(A,Xlefttemp)*Snew, Xlefttemp, Cold];
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
    end
    it_per_cycle(k) = j;    % final number of iterations for this cycle
    
	% Update solution
	Xleft = [Xleft, param.Vshort * W];
   
	% Compress solution
	[Xleft, S] = symm_compress_factors(Xleft, blkdiag(S, Snew),...
        param);
    solution_ranks(k) = size(Xleft, 2);

    min_core_eig(k) = min(eig(S));
    
	if res_approx < param.tol_res   % Clean the solution if it is indefinite
        flag = 1;
        if min(eig(S)) < 0
            [Xleft, S] = symm_compress_clean(Xleft, S);
            solution_ranks(k) = size(Xleft, 2);
        end
        break
    end
	
    % Set residual as next RHS
    C = [param.Vnext * param.Hnext,...
        param.Vshort * Y(:, end-r+1:end)];
    L = [zeros(r), eye(r); eye(r), zeros(r)];
    
    % Compress residual for restarting
    [C, L] = symm_compress_factors(C, L, param);

    r = size(C, 2);
    if param.verbose
        fprintf('i = %d, Rank of RHS: %d \n', k, r);
    end
    residual_ranks(k) = r;
end
if res_approx >= param.tol_res && k == param.max_restarts
    flag = 2;
end

% Update param outputs alphabetically
param.Acount = Acount;
param.flag = flag;
param.it_per_cycle = it_per_cycle;
param.matvecs = matvecs;
param.min_core_eig = min_core_eig;
param.num_restarts = k;
param.residual_norms = residual_norms;
param.residual_ranks = residual_ranks;
param.S = S;
param.solution_ranks = solution_ranks;
param.Xleft = Xleft;
end