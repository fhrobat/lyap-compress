function param = bgmres(A, B, param)
% param = BGMRES(A, B, param) computes the block GMRES approximation to
% the linear system AX = B.
%
% Possible configurations for param are described in the help line of
% BLOCK_SOLVE.  See also BLOCK_ARNOLDI.

%%
if param.verbose
    fprintf('%s_%s\n',mfilename,param.inner_product);
end

% Global variables and defaults
global Acount matvecs max_iterations
if isempty(Acount)
    Acount = 0;
end
if isempty(matvecs)
    matvecs = 0;
end
if isfield(param,'max_iterations')
    m = param.max_iterations;                                               % for independent runs
else
    m = max_iterations;                                                     % for use within EXTENDED_KRYLOV
end

s = size(B,2);
res_norm = 1;                                                               % first approximation is 0

if ~isfield(param,'res_scale')
    res_scale = block_norm(B, param.norm, A);                                % value by which to scale residual
else
    res_scale = param.res_scale;
end

for j = 1:m
    param = block_arnoldi(A,B,param);                                       % build basis, one step at a time
    Acount = Acount + 1;
    matvecs = matvecs + s;
    M = (param.Hshort'\param.Em)*param.Hnext'*param.Hnext;                  % compute GMRES modification matrix
    Xi = ((param.Hshort + M*param.Em')\param.E1)*param.H0;                  % compute solution to projected system

    res_norm(end+1) = block_norm([M; -param.Hnext]*...
            Xi(end-param.svec(end-1)+1:end,:),...
            param.norm, A)/res_scale;
    
    if param.verbose
        fprintf('          iter %d: residual norm = %g\n',...
            j,res_norm(end));
    end
    
    if res_norm(end) < param.tol_res
        param.Xm = param.Vshort*Xi;                                         % compute approximate solution
        param.res_norm = res_norm;
        param.Acount = Acount;
        param.matvecs = matvecs;
        return
    end
end
param.Xm = param.Vshort*Xi;                                                 % store current approximation

% Restart if no convergence and memory allowance exhausted
for k = 2:param.max_restarts
    param.M = M;
    B = [param.Vshort param.Vnext]*[param.M; -param.Hnext];                 % compute restart vector
    Gm = Xi(end-s+1:end,:);                                                 % cospatial factor
    
    param = clear_basis(param);
    
    for j = 1:m
        param = block_arnoldi(A, B, param);
        Acount = Acount + 1;
        matvecs = matvecs + s;
        
        M = (param.Hshort'\param.Em)*param.Hnext'*param.Hnext;              % compute GMRES modification matrix
        Xi = ((param.Hshort + M*param.Em')\param.E1)*param.H0*Gm;           % compute solution to projected system

        res_norm(end+1) = block_norm([M; -param.Hnext]*...
                Xi(end-param.svec(end-1)+1:end,:),...
                param.norm, A)/res_scale;
        
        if param.verbose
            fprintf('          iter %d: residual norm = %g\n',...
                j,res_norm(end));
        end
        
        if res_norm(end) < param.tol_res
            err_apx = param.Vshort*Xi;                                      % compute approximate solution
            param.Xm = param.Xm + err_apx;                                  % update solution
            param.Acount = Acount;
            param.matvecs = matvecs;
            return
        end
    end
end
    
param.Acount = Acount;
param.matvecs = matvecs;
param.res_norm = res_norm;
end