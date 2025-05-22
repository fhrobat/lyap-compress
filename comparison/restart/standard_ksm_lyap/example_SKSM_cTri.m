% Example: Solution of the Lyapunov equation 
%
% A*X+X*A+B*B^T=0,    A sym neg.def., B tall
%
% by the (standard) block Krylov subspace method (SKSM) with Galerkin condition,
% with enhanced computational strategy for the residual norm computation.
%
% Reference manuscript:
%
% Davide Palitta and Valeria Simoncini
% Computationally enhanced projection methods for symmetric Sylvester and Lyapunov matrix equations
% February 2016, ArXiv: 1602.05033
%
% Test example:
% A : 3600 x 3600 FD discretization of the Laplacian on [0,1]^2 w/ zero Dirichlet b.c.
% B : 3600 x 5 matrix w/ random entries (normal distribution), with unit Frobenius norm
% 
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS 
% FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
% COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER 
% IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 

nh=60;
T=diag(-2*ones(nh,1))+diag(ones(nh-1,1),1)+diag(ones(nh-1,1),-1);
I=speye(nh);
A=kron(T,I)+kron(I,T);
n=nh^2;
B=randn(n,5);
B=B/norm(B,'fro');
m=500;
tol=1e-6;


[Z,time_tot,time_res]=SKSM_Symm_Lyap(A,B,m,tol);
fprintf('\ntime_res: %10.5e\ntime_tot: %10.5e', time_res, time_tot)

% This matrix should never be formed for n large 
fprintf('\nFinal true absolute residual norm: %10.5e\n',norm((A*Z)*Z'+Z*(Z'*A')+B*B','fro'))
