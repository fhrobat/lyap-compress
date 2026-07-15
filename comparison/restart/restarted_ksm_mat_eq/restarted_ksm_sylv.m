function param = restarted_ksm_sylv(A, B, C, D, param)
% param = RESTARTED_KSM_SYLV(A, B, C, D, param) approximates the solution
% of Sylvester equation
%
%                       A * X + X * B + C * D' = 0
%
% via the restarted polynomial Krylov method with the parameters specified
% in the param struct.
%
% A detailed description of the variables and parameter specifications can
% be found in RESTARTED_KSM_MAT_EQ.
% 
% See also RESTARTED_KSM_LYAP and BLOCK_ARNOLDI.

%%
% Initialize
Acount = [0; 0];
flag = 0;
matvecs = [0; 0];
residual_ranks = [];
solution_ranks = [];
residual_norms = [];
it_per_cycle = [];

res_init = norm(C'*D,param.norm);
switch param.norm
        case 'fro'
            res_init = sqrt(trace((C'*C)*(D'*D)));   
        case {2, '2'}
            [~,RC]=qr(C,0);
            [~,RD]=qr(D,0);
            res_init = norm(RC*RD');
end    

r = size(C, 2);                                                             % assume same as D
Xleft = []; 
Xright = [];

% Set up distinct param structs for each space.
paramA = param;
paramB = param;

for k = 1:param.max_restarts
    if param.verbose == 1 || param.verbose == 2
        fprintf('Restart cycle: %d----------------------\n',k);
    end
    
    if k > 1    % Clear previous Krylov basis
        paramA = clear_basis(paramA);
        paramB = clear_basis(paramB);
    end
    
    if r == 0
        error('Compression tolerance is too high. Lower and try again.')
    end
    
    % Adjust memory allowance; assume each space has the same memory
    % limitations
    max_iterations = floor(param.memory_max/2/r);
    paramA.max_iterations = max_iterations;
    paramB.max_iterations = max_iterations;
    if max_iterations == 0
        error('Ran out of memory. Increase and try again.')
    end
    
	for j = 1:max_iterations
		% Augment the Krylov spaces
        paramA = block_arnoldi(A, C, paramA);
        paramB = block_arnoldi(B', D, paramB);
        Acount = Acount + [1; 1];
        matvecs = matvecs + r*[1; 1];

		% Solve the projected matrix equation
		Y = lyap(paramA.Hshort, paramB.Hshort',...
                paramA.E1 * paramA.H0 * (paramB.E1 * paramB.H0)');
        
        % Compute approximate residual
        switch param.norm
            case 'fro'
                res_approx =...
                    sqrt(norm(paramA.Hnext * Y(end-r+1:end,:),'fro')^2 ...
                    + norm(Y(:,end-r+1:end) * paramB.Hnext','fro')^2)/...
                    res_init;
            case {2, '2'}
                res_approx =...
                    max(norm(paramA.Hnext * Y(end-r+1:end,:), 2),...
                    norm(Y(:,end-r+1:end) * paramB.Hnext', 2))/...
                    res_init;
        end
		
        % Print residual
        switch param.verbose
            case 0    % Don't print anything
                residual_norms = [residual_norms, res_approx];
                
            case 1  % Print the approximate residual
                fprintf('     Inner iter: %d, res = %e\n', j, res_approx);
                residual_norms = [residual_norms, res_approx];
                
            case 2   % Calculate and print the true residual
                X = paramA.Vshort * Y * paramB.Vshort';
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
        
    end
    it_per_cycle(k) = j;    % final number of iterations for this cycle

	% Update the solution
	[W, S, Z] = svd(Y);
	Xleft = [Xleft, paramA.Vshort * W * sqrt(S)];
	Xright = [Xright, paramB.Vshort * Z * sqrt(S)];
    
	% Compress the solution
	[Xleft, Xright] = compress_factors(Xleft, Xright, param);
	solution_ranks(k) = size(Xleft, 2);

	if res_approx < param.tol_res
        flag = 1;
		break
	end
	
	% Set residual as next RHS
	C = [paramA.Vnext * paramA.Hnext,...
        paramA.Vshort * Y(:,end-r+1:end)];
	D = [paramB.Vshort * Y(end-r+1:end,:)',...
        paramB.Vnext * paramB.Hnext];
	
	% Compress residual for restarting
	[C, D] = compress_factors(C, D, param);

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
param.num_restarts = k;
param.residual_norms = residual_norms;
param.residual_ranks = residual_ranks;
param.solution_ranks = solution_ranks;
param.Xleft = Xleft;
param.Xright = Xright;
end