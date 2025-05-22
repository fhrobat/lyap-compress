function param = extended_ksm_mat_eq(A,B,C,D,param)
% param = EXTENDED_KSM_MAT_EQ(A,B,C,D,param) is a wrapper function for
% calling extended Krylov subspace methods for the matrix equation
%
%                       A * X + X * B + C * D' = 0.
%
% Possible configurations are as follows:
%
% param = EXTENDED_KSM_MAT_EQ(A,[],C,[],param) calls EXTENDED_LYAP for the
% Lyapunov problem.
%
% param = EXTENDED_KSM_MAT_EQ(A,[],C,D,param)
%   and
% param = EXTENDED_KSM_MAT_EQ(A,B,C,[],param)
%   and
% param = EXTENDED_KSM_MAT_EQ(A,B,C,D,param) call EXTENDED_SYLV for the
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
% A and B do not have to be numeric. See IMPLICIT_STRUCT for how to
% construct operators with implicit multiplication and inversion
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
%                           unless tolerance satisfied
%
%   .memory_max             maximum number of vectors allowed in the memory
%                           buffer; used to dynamically determine bases
%                           size
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
%   .residual_norms     vector of norms of the residual matrix per
%                       iteration
%
%   .S                  core matrix such that Xleft*S*Xright'
%                       (Xleft*S*Xleft' for symmetric Lyapunov) is the
%                       computed solution
%
%   .Xleft              low-rank factorization of the solution; only Xleft
%                       is generated for the fully symmetric Lyapunov
%                       problem to give X = Xleft*S*Xleft'
%
%   .Xright             low-rank factorization of the solution; X =
%                       Xleft*Xright'
%
% See also RESTARTED_KSM_MAT_EQ, STANDARD_KSM_LYAP, EXTENDED_KSM_LYAP, and
% EXTENDED_KSM_SYLV.

%%
addpath(genpath('extended_ksm_mat_eq'))
if isnumeric(A)
    % Turn A into a struct with precomputed Cholesky or LU factors
    A = struct_make(A);
end
if exist('B','var') && ~isempty(B)
    if isnumeric(B)
        B = struct_make(B);
    end
end
    

% Fill-in or initialize parameter struct
if ~exist('param','var') || isempty(param)
    param = param_init_mat_eq;
else
    param = param_init_mat_eq(param);
end

if isempty(B) && isempty(D)         % Lyapunov
    param = extended_ksm_lyap(A,C,param);
else                                % Sylvester
    param = extended_ksm_sylv(A,B,C,D,param);
end
end