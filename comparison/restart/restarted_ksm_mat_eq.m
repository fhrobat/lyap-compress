function param = restarted_ksm_mat_eq(A,B,C,D,param)
% param = RESTARTED_KSM_MAT_EQ(A,B,C,D,param) is a wrapper function for
% calling the restarted Lyapunov and restarted Sylvester routines to solve
%
%                       A * X + X * B + C * D' = 0.
%
% Possible configurations are as follows:
%
% param = RESTARTED_KSM_MAT_EQ(A,[],C,[],param) calls RESTARTED_LYAP for
% the Lyapunov problem.
%
% param = RESTARTED_KSM_MAT_EQ(A,[],C,D,param)
%   and
% param = RESTARTED_KSM_MAT_EQ(A,B,C,[],param)
%   and
% param = RESTARTED_KSM_MAT_EQ(A,B,C,D,param) call RESTARTED_SYLV for the
% Sylvester problem.
%
%-------------- INPUTS ----------------------------------------------------
% A             linear coefficient of the equation; the eigenvalues of A
%               must have negative real part
%
% B             (optional) second linear coefficient of the Sylvester
%               equation; negative eigenvalues must be disjoint from the
%               spectrum of A
%
% A and B do not have to be numeric. See STRUCT_MAKE and PRECOND_STRUCT for
% how to construct operators with implicit multiplication and inversion
% procedures.
%
% C             low-rank factorization of the right-hand side
%
% D             low-rank factorization of the right-hand side for the
%               nonsymmetric Lyapunov or Sylvester problem
%
% param         a struct containing parameter specifications:
%
%   .max_restarts           maximum number of restarts to be performed
%                           unless tolerance satisfied; only applies to
%                           outer iteration
%
%   .memory_max             maximum number of vectors allowed in the memory
%                           buffer; used to dynamically determine bases
%                           sizes
%
%   .norm                   choice of norm for measuring residual; possible
%                           options are the same as for built-in norm
%
%   .tol_res                residual tolerance for the stopping criterion
%
%   .tol_comp               compression tolerance
%
%   .verbose                if == 0, nothing prints; if == 1, prints the
%                           approximate residual at each step; if == 2,
%                           prints the approximate and true residual
%
%-------------- OUTPUTS ---------------------------------------------------
% param         the same struct updated to contain outputs
%
%   .Xleft              low-rank factorization of the solution; only Xleft
%                       is generated for the fully symmetric Lyapunov
%                       problem to give X = Xleft*Xleft'
%
%   .Xright             low-rank factorization of the solution; X =
%                       Xleft*Xright'
%
%   .S                  core matrix such that Xleft*S*Xright'
%                       (Xleft*S*Xleft' for symmetric Lyapunov) is the
%                       computed solution
%
%   .residual_ranks     vector of ranks of the residual matrix per
%                       restart cycle
%
%   .solution_ranks     vector of ranks of the approximate solution per
%                       restart cycle
%
%   .residual_norms     vector of norms of the residual matrix per
%                       iteration
%
%   .num_restarts       number of restarts performed
%
%   .Acount             a vector of the number of times the matrix A is
%                       accessed per iteration; when B is also specified,
%                       the second row corresponds to the same vector but
%                       for B
%
%   .flag               a double 0, 1, or 2, indicating whether the method
%                       has converged
%
%   .matvecs            a vector of the number of vectors or columns at
%                       which A (respectively, B) is evaluated
%
% See also EXTENDED_KSM_MAT_EQ, STANDARD_KSM_LYAP, RESTARTED_KSM_LYAP, and
% RESTARTED_KSM_SYLV.

%%
addpath(genpath('restarted_ksm_mat_eq'));

% Fill-in or initialize parameter struct
if ~exist('param','var') || isempty(param)
    param = param_init_mat_eq;
else
    param = param_init_mat_eq(param);
end

if isempty(B) && isempty(D)         % Lyapunov
    param = restarted_ksm_lyap(A,C,param);
else                                % Sylvester
    param = restarted_ksm_sylv(A,B,C,D,param);
end
