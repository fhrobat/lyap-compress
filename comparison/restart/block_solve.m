function [X, param] = block_solve(A,B,param)
% param = BLOCK_SOLVE(A,B,param) is a wrapper function for calling block
% Krylov subspace routines to solve A*X = B.
%
% Possible configurations are as follows:
%   X = block_solve(A,B) solves A*X = B for default parameter settings
%   determined by PARAM_INIT_BLOCK_SOLVE.
%
%   X = block_solve(A,B,param) solves A*X = B for the provided parameter
%   settings.
%
%   [X, param] = block_solve(A,B,param) solves A*X = B for the provided
%   parameter settings and returns the same param struct, now packed with
%   outputs and updated parameter settings.
%
%-------------- INPUTS ----------------------------------------------------
% A             matrix or struct with a .multiply field defined as the
%               anonymous function @(rho, eta, x) rho * A * x - eta * x. If
%               A is a struct and non-Hermitian, then it must also have a
%               .trans.multiply field defined as the anonymous function
%               @(rho, eta, x) rho * A * x - eta * x.
%
% B             right-hand side vector or block vector (i.e.,
%               tall-and-skinny matrix)
%
% param         a struct containing parameter specifications:
%
%   .bcg_variant            a string specifying which retooled version of
%                           BCG to use; options include 'standard', 'bcga',
%                           'bcgadq', or 'bcgrq'; default is 
%                           'standard'. For a description of the variants,
%                           see BCG.
%
%   .max_iterations         maximum number of iterations allowed; note that
%                           one block basis vector of size(B) is generated
%                           per iteration
%
%   .max_restarts           maximum number of restarts allowed (for BGMRES;
%                           does not apply to BCG)
%
%   .norm                   a string specifying which norm to use for
%                           the error measure; default is 'fro'
%
%   .precond_left           a left preconditioner matrix L (i.e., solve
%                           instead L\(A*X) = L\B)
%
%   .precond_right          a right preconditioner matrix R (i.e., solve
%                           instead A*(R\Y) = B, X = R\B)
%                           If both preconditioners are specified, both
%                           will be used. The default is no
%                           preconditioning.
%
%   .solver                 a string specifying which solver to use;
%                           options include 'bcg', 'bfom', and 'bgmres'.
%                           The default for hermitian matrices is 'bcg';
%                           otherwise 'bgmres'.
%
%   .tol_res                residual tolerance
%
%   .verbose                if == 0, nothing prints; if == 1, prints the
%                           specified convergence check measure at each
%                           step
%
%-------------- OUTPUTS (FINAL) -------------------------------------------
% X                         the approximation solution
%
% param         the same struct updated to contain outputs
%
%   .res_norm               a vector of residual norms measured per
%                           iteration if .inner_cycle = true; otherwise per
%                           restart cycle
%
%-------------- OUTPUTS (INTERMEDIATE) ------------------------------------
% param         the same struct updated to contain outputs
%
%   .H0                     the block norm of the vector used to build the
%                           Krylov basis
%
%   .Hshort                 square upper Hessenberg matrix satisfying
%                           Arnoldi/Lanczos relation
%
%   .Hnext                  small block matrix completing .Hshort
%
%   .Vshort                 a concatenation of m basis vectors
%
%   .Vnext                  the m+1 basis vector
%
%   .E1                     block unit vector with I in the first block
%
%   .Em                     block unit vector with I in the last block
%
%   .M                      modification matrix (only in BGMRES)
%
% See also BCG, BGMRES, STRUCT_MAKE, or PRECOND_STRUCT.

%%
addpath(genpath('solvers'))
param = param_init_block_solve(param);                                      % set defaults for param

% To handle preconditioning, we turn numeric A into a struct with just a
% multiply field. If A is a struct has been built by STRUCT_MAKE, then
% PRECOND_STRUCT will be accessed whenever A.solve is called.
if isnumeric(A) &&...
        (isfield(param,'precond_left') || isfield(param,'precond_right'))
    A = precond_struct(A, param);
end

param.solver = lower(param.solver);
switch param.solver
    case {'cg', 'bcg'}
        if isfield(param,'precond_left')
            param = bcg(A, param.precond_left\B, param);
        else
            param = bcg(A, B, param);
        end
        
    case {'gmres', 'bgmres'}
        if isfield(param,'precond_left')
            param = bgmres(A, param.precond_left\B, param);
        else
            param = bgmres(A, B, param);
        end
end
% Recover solution from right preconditioned problems
if isfield(param,'precond_right')
    param.Xm = param.precond_right\param.Xm;
end

X = param.Xm;
param = rmfield(param, 'Xm');

end