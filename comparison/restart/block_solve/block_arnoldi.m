function param = block_arnoldi(A,B,param)
% param = BLOCK_ARNOLDI(A,B,param) computes the block Arnoldi basis with
% respect to A, B, and the classical inner product without deflation. One
% must assume that exact deflation will not occur in order to use this
% algorithm; there are no fail-safes for breakdowns.
%
% Possible configurations for param are described in the help line of
% BLOCK_SOLVE.  See also BGMRES.

%%
[n, s] = size(B);

if ~isfield(param,'Hshort')
    % Start the basis from scratch
    param.Acount = 0;
    param.matvecs = 0;
    
    jnew = 1:s;
    V = zeros(n, 2*s);
    [V(:,1:s), param.H0] = qr(B,0);
    H = zeros(2*s,s);
    
    W = struct_mult(A,V(:,jnew));
    param.Acount = param.Acount + 1;
    param.matvecs = param.matvecs + s;
    
    jold = 1:s;
    H(jold,jnew) = V(:,jold)'*W;
    W = W - V(:,jold)*H(jold,jnew);
    jold = jnew;
    jnew = jold(end)+1:jold(end)+s;
    [V(:,jnew), H(jnew,jold)] = qr(W,0);    
    
    param.svec = s*ones(1,2);
    
else
    j = size(param.Hshort,1)/s;
    
    % Continue the basis
    V = zeros(n, (j+2)*s);
    V(:,1:j*s) = param.Vshort;
    V(:,j*s+1:(j+1)*s) = param.Vnext;
    H = zeros((j+2)*s,(j+1)*s);
    H(1:j*s,1:j*s) = param.Hshort;
    H(j*s+1:(j+1)*s, (j-1)*s+1:j*s) = param.Hnext;

    j = j+1;
    jnew = (j-1)*s+1:j*s;

    W = struct_mult(A,V(:,jnew));
    param.Acount = param.Acount + 1;
    param.matvecs = param.matvecs + s;
    for jj = 1:j
        jold = (jj-1)*s+1:jj*s;
        H(jold,jnew) = V(:,jold)'*W;
        W = W - V(:,jold)*H(jold,jnew);
    end
    jold = jnew;
    jnew = jold(end)+1:jold(end)+s;
    [V(:,end-s+1:end) , H(jnew,jold)] = qr(W,0);
    
    param.svec = s*ones(1,j+1);
end

Is = eye(s);
E1 = zeros(size(H,2),s); E1(1:s,1:s) = Is;
param.E1 = E1;

Em = zeros(size(H,2),s);
Em(end-s+1:end,end-s+1:end) = Is;
param.Em = Em;

param.Hnext = H(end-s+1:end,end-s+1:end);
param.Hshort = H(1:end-s,:);
param.Vnext = V(:,end-s+1:end);
param.Vshort = V(:,1:end-s);
end