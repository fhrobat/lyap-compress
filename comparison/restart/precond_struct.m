function S = precond_struct(A, param, trans)
% S = PRECOND_STRUCT(A, param) returns a struct with a multiply field
% and a preconditioned version of the matrix or struct A.  If A is a
% struct, it must have a .multiply field defined as the anonymous function
%
%                @(x) A * x.
%
% S = PRECOND_STRUCT(A, param, trans), with trans == true, returns a
% struct with just a multiply field and a preconditioned version of A'.  If
% A is a struct, then it must have a .trans.multiply field defined as the
% anonymous function
%
%                @(x) A' * x.
%
% See also STRUCT_MAKE and BLOCK_SOLVE.

%%
if nargin == 2
    trans = false;
end

% Set up local matrices for preconditioners.  Preconditioners are set to
% identity if not specified by param.
if isfield(param,'precond_left') && isfield(param,'precond_right')
    L = param.precond_left;
    R = param.precond_right;
elseif isfield(param,'precond_left') && ~isfield(param,'precond_right')
    L = param.precond_left;
    R = speye(size(A));
elseif ~isfield(param,'precond_left') && isfield(param,'precond_right')
    L = speye(size(A));
    R = param.precond_right;
elseif ~isfield(param,'precond_left') && ~isfield(param,'precond_right')
    L = speye(size(A));
    R = speye(size(A));
end

% Note that we have to be pedantic with parentheses, otherwise MATLAB
% computes blindly from left to right, which is undesirable.
if ~trans
    if isnumeric(A)
        S = struct('multiply', @(x) L\ (A* (R\x) ) );
    elseif isstruct(A)
        S = struct('multiply', @(x) L\ A.multiply(R\x) );
    end
elseif trans
    if isnumeric(A)
        S = struct('multiply', @(x) R'\ (A'* (L'\x) ) );
    elseif isstruct(A)
        S = struct('multiply', @(x) R'\ A.trans.multiply(L'\x) );
    end
end
end
