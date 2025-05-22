function W = struct_mult(A,V)
% W = STRUCT_MULT(A,V) returns the action of A on the matrix V.  A must be
% be either a matrix or a struct formatted as in STRUCT_MAKE.
%
% See also PRECOND_STRUCT.

%%
if isnumeric(A)
    W = A*V;
elseif isstruct(A)
    W = A.multiply(V);
end
end