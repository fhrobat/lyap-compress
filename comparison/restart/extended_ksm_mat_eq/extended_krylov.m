function param = extended_krylov(A, B, param)
% param = EXTENDED_KRYLOV(A, B, param) builds an extended Krylov basis
% given the operator A and starting (block) vector B, according to the
% specifications in struct param. The number of A-calls (Acount) and the
% number of vectors to which A is applied (matvecs) are stored as
% param.Acount and param.matvecs, respectively, and cleared at the
% beginning of each run.
%
% Notes:
% --When A is numeric, backslash is used for inversion and is
% counted as a single operation with A.
% --Basis vectors are added in the order inv(A)*B then A*B, and they are
% added two at a time.
%
% The basis vectors are stored as param.V, and projection matrices as
% param.H and param.K, where H and K are such that
%
%                 A*V*K = V*H
%
% For further param options, see BLOCK_SOLVE.  For how to build an operator
% with implicit solves, see STRUCT_MAKE and PRECOND_STRUCT.

%%
global Acount matvecs max_iterations
if ~strcmp(A.solver,'backslash')
    max_iterations = param.max_iterations;
end
Acount = 0;
matvecs = 0;
s = size(B,2);

if ~isfield(param,'H')
    % Start the basis from scratch
    [V, param.H0] = qr(B, 0);
    
    H = zeros(s, 0);
    K = zeros(s, 0);
    
    [V, K, H, W] = add_zero_pole(V, K, H, A, V);
    [param.V, param.K, param.H, param.last] = add_inf_pole(V, K, H, A, W);
    
else
    % Continue the basis
    [V, K, H, W] =...
        add_zero_pole(param.V, param.K, param.H, A, param.last);
    [param.V, param.K, param.H, param.last] =...
        add_inf_pole(V, K, H, A, W);
end
param.Acount = Acount;
param.matvecs = matvecs;

Is = eye(s);
E1 = zeros(size(param.H,2),s); E1(1:s,1:s) = Is;
param.E1 = E1;
Em = zeros(size(param.H,2),s);
Em(end-s+1:end,end-s+1:end) = Is;
param.Em = Em;
end

%% Subfunctions
function [V, K, H, W] = add_zero_pole(V, K, H, A, W)
% ADD_ZERO_POLE(V, K, H, A, W) is an utility routine that adds a zero pole
% to the space. The vector W is the continuation vector.
    global Acount  matvecs
    s = size(W, 2);
    
    W = A.solve(W);
    if strcmp(A.solver, 'backslash')
        Acount = Acount + 1;
        matvecs = matvecs + s;
        
        % Otherwise, Acount and matvecs calculated globally in bgmres or
        % bcg
    end
    H0 = V'*W;
    W = W - V*H0;

    % Enlarge H and K
    H(size(H, 1) + s, size(H, 2) + s) = 0;
    K(size(K, 1) + s, size(K, 2) + s) = 0;

    K(1:end-s, end-s+1:end) = H0;
    H(end-2*s+1:end-s, end-s+1:end) = eye(s);

    [W, R] = qr(W, 0);    
    K(end-s+1:end, end-s+1:end) = R;

    V = [V, W];
end

function [V, K, H, W] = add_inf_pole(V, K, H, A, W)
% ADD_INF_POLE is a utility routine that adds an infinity pole to the
% space. The vector W is the continuation vector.
    global Acount  matvecs
    s = size(W, 2);
    
    W = A.multiply(W);
    Acount = Acount + 1;
    matvecs = matvecs + s;
    H0 = V'*W;
    W = W - V*H0;

    % Enlarge H and K
    H(size(H, 1) + s, size(H, 2) + s) = 0;
    K(size(K, 1) + s, size(K, 2) + s) = 0;

    H(1:end-s, end-s+1:end) = H0;
    K(end-2*s+1:end-s, end-s+1:end) = eye(s);

    [W, R] = qr(W, 0);    
    H(end-s+1:end, end-s+1:end) = R;

    V = [V, W];
end

