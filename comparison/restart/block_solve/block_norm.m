function N = block_norm(X,normstr,A)
% N = BLOCK_NORM(X,normstr,A) computes a scalar norm for the block vector X

%%
switch normstr
    case 'A-fro'                                                            % A-weighted Frobenius norm
        N = sqrt(real(trace(X'*(struct_mult(A,X)) ) ) );
    case 'fro'                                                              % Frobenius norm; default
        N = norm(X,'fro');
    case 'first-col'
        N = norm(X(:,1));
    case 'last-col'
        N = norm(X(:,end));
    case {2,'2'}                                                            % matrix 2-norm
        N = norm(X,2);
    case {1,'1'}                                                            % matrix 1-norm
        N = norm(X,1);
    case {inf,Inf,'inf','Inf'}                                              % matrix inf-norm
        N = norm(X,Inf);
end