function param = standard_ksm_lyap(A, C, param)
% param = STANDARD_KSM_LYAP(A, C, param) returns the solution of the
% Lyapunov equation
%
%                    A * X + X * A' + C * C' = 0
%
% where A is symmetric negative definite, and C is tall and skinny, by the
% (standard) block Krylov subspace method (SKSM) with Galerkin condition
% and enhanced computational strategy for the residual norm computation.
% 
%-------------- INPUTS ----------------------------------------------------
% A             symmetric negative definite matrix
%
% A does not have to be numeric. See PRECOND_STRUCT for how to construct
% operators with implicit multiplication and inversion procedures.
%
% C             right-hand side matrix
%
% param         a struct containing parameter specifications:
%
%   .max_iterations         maximum number of basis vectors to generate;
%                           note that we have no memory restrictions here,
%                           thanks to short-term recurrences
%
%   .tol_res                residual tolerance for the stopping criterion
%                           in Frobenius norm
%
%   .verbose                if == 0, nothing prints; if == 1, prints the
%                           approximate residual at each step
%
%-------------- OUTPUTS ---------------------------------------------------
% param         the same struct updated to contain outputs
%
%   .Acount             a vector of the number of times the matrix A is
%                       accessed per iteration
%
%   .flag               a double 0, 1, or 2, indicating whether the method
%                       has converged
%
%   .matvecs            a vector of the number of vectors or columns at
%                       which is evaluated
%
%   .residual_norms     vector of norms of the residual matrix per
%                       iteration
%
%   .Xleft              tall-and-skinny matrix such that Xleft*Xleft' is
%                       the approximate solution (Xleft is truncated with
%                       compression tolerance 1e-12)
%
% See also RESTARTED_KSM_MAT_EQ and EXTENDED_KSM_MAT_EQ.
%
% -------------------------------------------------------------------------
% Reference manuscript:
%
% 	Davide Palitta and Valeria Simoncini. Computationally enhanced
% 	projection methods for symmetric Sylvester and Lyapunov matrix
% 	equations. Journal of Computation and Applied Mathematics, Vol. 330,
% 	pg. 648-659, 2018.
%
% email: davide.palitta3@unibo.it, valeria.simoncini@unibo.it
%
% REQUIREMENT: the subroutines EIG_DSTEVR and TRIDIAG_DSBTRD which use
% LAPACK and C-BLAS subroutines.
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.
%
% See https://zenodo.org/record/3252320#.XjFKGchKg2x for the original SKSM
% routines.
%
% Modified November 2019 by Kathryn Lund.

%%
addpath(genpath('standard_ksm_lyap'))
if nargin == 2
    param = param_init_mat_eq;
    param.max_iterations = 100;
end

% Initialize
Acount = 0;
flag = 0;
matvecs = 0;
residual_norms = [];
m = param.max_iterations;

res_init = norm(C'*C,param.norm);                                           % scale by initial residual

r = size(C,2);

% Orthonormalize C
[U, beta] = qr(C,0);
beta2 = beta*beta';

% Preallocation
TA = zeros((m+1)*r,m*r);
TTA = zeros(r+1,m*r);

for j = 1:m
    jms = (j-1)*r+1;
    j1s = (j+1)*r;
    js = j*r;
    js1 = js+1;
   
    % Sequence in A
    % New basis block (double modified Gramm-Schmidt)
    Up = struct_mult(A, U(:, end-r+1:end));
    Acount = Acount + 1;
    matvecs = matvecs + r;
  
    for l = 1:2
        K = min([j, 3]);
        for kk = 1:K
              k1 = (kk-1)*r+1;
              k2 = kk*r;
              coeff = U(:, k1:k2)' * Up;
              TA((j-K+kk-1)*r+1:(j-K+kk)*r, jms:js) =...
                  TA((j-K+kk-1)*r+1:(j-K+kk)*r, jms:js) + coeff; 
              Up = Up - U(:, k1:k2) * coeff; 
        end
    end
  
    if j == 1
        [U(:, r+1:2*r), TA(js1:j1s, jms:js)] = qr(Up, 0);
        
        % Prepare LLA for dstrd
        for i = 1:r
            TTA(end-i+1, i:r) = diag(TA(jms:js, jms:js), i-1);
        end
        LLA = TTA(:, 1:j*r);
        
    else
        % Move the blocks
        if j > 2
            U(:, 1:2*r) = U(:, r+1:3*r);
        end
        [U(:, 2*r+1:3*r), TA(js1:j1s, jms:js)] = qr(Up,0);
    
        % Prepare LLA for dstrd
        HHA = zeros(r+1,r);
        for i = 1:r+1
            HHA(end-i+1,1:r) = diag(TA(js-r+2-i:js-i+1, js-r+1:js));
        end
        TTA(:,(j-1)*r+1:j*r) = HHA;
        LLA = TTA(:,1:j*r);
    
    end

    % Take the ((j+1)s,js) block of TA 
    ll = TA(js1:j1s, jms:js)';

    % Tridiagonalize LLA by dsbtrd
    [D, E, P] = tridiag_dsbtrd(LLA);
    
    % Compute the eigendecompositions of the obtained tridiagonal matrices
    % by dstevr
    [Q, lambda] = eig_dstevr(D,E);
    
    
    % Compute the residual norm
    S = ((P(1:r,:) * Q)' * beta2 * (P(1:r,:)*Q))./...
        (lambda * ones(1,js) + ones(js,1) * lambda');
    W = (P(end-r+1:end, :) * Q)' * ll;
    
    res_approx = sqrt(2)*norm(S*W,'fro')/res_init;
    residual_norms(end+1) = res_approx;
    
    if param.verbose
        fprintf('Iter: %d, res = %e\n', j, res_approx)
    end
    
    if res_approx < param.tol_res
        flag = 1;
        break
    end
end
it = j;

% Compute the full eigendecomposition of T and the matrix Y
TA(1:js,1:js) = (TA(1:js, 1:js) + TA(1:js, 1:js)') / 2;
[QA, lambda] = eig(full(TA(1:js, 1:js)));
lambdaA = diag(lambda);
Q = QA(1:r,:);

% Compute the solution of the last projected equation
Y = -(Q'*beta2*Q) ./ (lambdaA*ones(1,js) + ones(js,1)*lambdaA');

% Truncation strategy
[uY, sY] = eig(Y);
[sY, id] = sort(abs(diag(sY)));
sY = flipud(sY);
uY = uY(:,id(end:-1:1));
is = sum(abs(sY)>1e-12);
Y0 = QA*(uY(:,1:is) * diag(sqrt(sY(1:is))));

% Two-pass
[U, ~] = qr(C,0);
Z = U * Y0(1:r,:);

% separate j=1,2 from the other j's to avoid "if" within the "for" loop
j = 1;
j1s = 2*r;
js1 = r + 1;

% sequence in A
UpA = struct_mult(A, U(:, end-r+1:end));
Acount = Acount + 1;
matvecs = matvecs + r;
U(:, r+1:2*r) = (UpA - U(:, end-r+1:end)*...
    TA((j-1)*r+1:j*r, (j-1)*r+1:j*r))/...
    TA(j*r+1:(j+1)*r, (j-1)*r+1:j*r);
Z = Z + U(:,end-r+1:end) * Y0(js1:j1s,:);

j = 2;
j1s = 3*r;
js1 = 2*r+1;

% sequence in A
UpA = struct_mult(A, U(:,end-r+1:end));
Acount = Acount + 1;
matvecs = matvecs + r;
U(:,2*r+1:3*r) = (UpA-U(:,1:2*r)*...
    TA((j-2)*r+1:j*r, (j-1)*r+1:j*r))/...
    TA(j*r+1:(j+1)*r, (j-1)*r+1:j*r);
Z = Z + U(:,end-r+1:end) * Y0(js1:j1s,:);

for j = 3:it-1
    j1s = (j+1)*r;
    js1 = j*r+1;
    
    % sequence in A
    UpA = struct_mult(A, U(:, end-r+1:end));
    Acount = Acount + 1;
    matvecs = matvecs + r;
    U(:, 1:2*r) = U(:, r+1:3*r);
    U(:, 2*r+1:3*r) = (UpA - U(:, 1:2*r)*...
        TA((j-2)*r+1:j*r, (j-1)*r+1:j*r))/...
        TA(j*r+1:(j+1)*r, (j-1)*r+1:j*r);
    Z = Z + U(:, end-r+1:end) * Y0(js1:j1s,:);
end

param.Acount = Acount;
param.flag = flag;
param.matvecs = matvecs;
param.residual_norms = residual_norms;
param.Xleft = Z;