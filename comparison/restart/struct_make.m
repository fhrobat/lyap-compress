function S = struct_make(A, param)
% S = STRUCT_MAKE(A, param) returns a struct with fields specifying how to
% multiply and solve with A, with specifications for a block Krylov
% subspace method provided by param; see BLOCK_SOLVE for details.The block
% Krylov solver should declare global variables Acount and matvecs to count
% the number of calls to the operator A and the number of columns to which
% A is applied, respectively.  Left and right preconditioners should be
% stored as param.precond_left and param.precond_right, respectively.
% BLOCK_SOLVE will automatically incorporate them.
%
% If only A is provided, then a Cholesky/LU factorization is computed and
% stored.  This default is labeled 'backslash' in the solver field, as it
% is what the built-in Matlab operator does.  However, by precomputing and
% storing factors, speed-ups are possible if A.solve is called repeatedly.
%
% The output S is a struct with the following fields:
%
%   .multiply(x)    an inline function taking a (block) vector x and
%                   returning A * x
%
%   .solve(x)       an inline function taking a (block) vector x and
%                   returning an approximation y to the system A y = x
%
%   .solver         a string specifying which solver is used, e.g., 'bcg',
%                   'cg', 'bgmres', 'gmres', or 'backslash'.
%
%   .trans          a struct with the same fields (.multiply, .solve,
%                   .solver) as the main struct, but corresponding to the
%                   action of the transpose of A
%
% See also STRUCT_MULT, PRECOND_STRUCT, and PARAM_INIT_BLOCK_SOLVE.

%%
if nargin == 1  % then precompute Cholesky/LU factors
    [R, flag, P] = chol(A, 'vector');             % R'*R = P'*A*P
    sgn = 1;
    
    % Check if negative definite
    if flag
        [R, flag, P] = chol(-A, 'vector');
        sgn = -1;
    end
    tmp = P(1:size(A,1));
    IP(tmp) = 1: size(A,1); % Inverse permutation
    if ~flag
        S = struct(...
	    'solve', @(x) sparse_chol_solve(R, P, IP, x, sgn) , ...
            'multiply', @(x) A * x,...
            'solver', 'backslash');
        S.trans = S;
        
    else
        [L, U, P, Q] = lu(A, 'vector');                 % L*U = P*A*Q  
        tmp = P(1:size(A, 1));
    	IP(tmp) = 1:size(A, 1);
        tmp = Q(1:size(A, 1));
    	IQ(tmp) = 1:size(A, 1);
        S = struct(...
            'solve', @(x) sparse_lu_solve(L, U, P, IQ, x),...
            'multiply', @(x) A * x,...
            'solver', 'backslash');

        S.trans = struct(...            % U'*L' = P' * A'
            'solve', @(x) sparse_trasp_lu_solve(L, U, IP, Q, x),...
            'multiply', @(x) A' * x,...
            'solver', 'backslash');
    end
    
else
    if ~isfield(param, 'hermitian')
        param.hermitian = false;
    end
    if ~isfield(param, 'symmetric')
        param.symmetric = false;
    end
    if param.hermitian || param.symmetric
        S = struct(...
            'solve', @(x) block_solve(A, x, param),...
            'multiply', @(x) A * x,...
            'solver', param.solver);
        
        S.trans = S;
        
    else
        S = struct(...
            'solve', @(x) block_solve(A, x, param),...
            'multiply', @(x) A * x,...
            'solver', param.solver);

        % Set up new param struct for the transpose
        trans_param = param;
        if isfield(param,'precond_left') && isfield(param,'precond_right')
            trans_param.precond_right = param.precond_left';
            trans_param.precond_left = param.precond_right';
        elseif isfield(param,'precond_left') && ~isfield(param,'precond_right')
            trans_param.precond_right = param.precond_left';
            trans_param = rmfield(trans_param,'precond_left');
        elseif ~isfield(param,'precond_left') && isfield(param,'precond_right')
            trans_param = rmfield(trans_param,'precond_right');
            trans_param.precond_left = param.precond_right';
        end

        S.trans = struct(...
            'solve', @(x) block_solve(A', x, trans_param),...
            'multiply', @(x) A' * x,...
            'solver', param.solver);
    end
    S.param = param;                                                        % necessary for updating parameters later
end
end

%% Auxiliary functions
function y = sparse_chol_solve(R, P, IP, x, sgn)
    y = sgn * R\(R'\x(P, :));
    y = y(IP, :);	
end
function y = sparse_lu_solve(L, U, P, IQ, x)
	y = U \ (L \ x(P, :));
	y = y(IQ, :);
end
function y = sparse_trasp_lu_solve(L, U, IP, Q, x)
	y = L' \ (U' \ x(Q, :));
	y = y(IP, :);
end
