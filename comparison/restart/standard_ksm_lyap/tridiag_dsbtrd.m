function [varargout] = tridiag_dsbtrd(varargin)
% TRIDIAG_LAPACK computes the matrix tridiagonal decomposition 
% by calling LAPACK subroutines.
%
% [Q, T] = tridiag_lapack (A)
% [Q, T] = tridiag_lapack (A,UPLO)
% R = tridiag_lapack (A, ...)
% [Q, D, E] = tridiag_lapack (A, ...)
% [Q, D, E, err] = tridiag_lapack (A, ...)
%
% Compute a real symmetric tridiagonal matrix decomposition of 
% a complex Hermitian matrix A:
%          A = Q*T*Q'
% where A is complex Hermitian or real symmetric, 
% Q is unitary or orthogonal, and
% T is real symmetric tridiagonal.
%
% If the second argument UPLO is used, the order of the reflector 
% vectors are chosen according to UPLO value 'L' (default) or 'U'.
% See http://www.netlib.org/lapack/explore-html/df/d7d/zhetrd_8f_source.html
% 
% The tridiagonalization has applications in matrix reduction 
% preserving eigenvalues. Thus, the eigenvalues of T are 
% identical to those of A. 
%
% When called with one return value, it returns R, a compressed form of T 
% and Householder reflector of Q, which is computed in LAPACK 
% subroutine DSBTRD.
% 
% When called with 3 return values, it returns D and E, which are
% the diagonal vector and the off-diagonal vector of T:
%          T = diag(D) + diag(E,1) + diag(E,-1)
%
% When called with 4 return values, it returns the residual error:
%          err = norm(A-Q*T*Q')/norm(A).
%
% lapack interface routines come from the work of Tim Toolan.
% See http://www.mathworks.co.kr/matlabcentral/fileexchange/16777-lapack

if nargin<1,
	varargout = {};
	help tridiag_dsbtrd;
	return;
end
if nargin>0,
	A = varargin{1};
end
if nargin>1,
	UPLO = varargin{2};
else
	UPLO = 'U';
end
if UPLO~='U' && UPLO~='L'
	varargout = {[]};
	disp('> UPLO must be ''L'' or ''U''.');
	return;
end


% lapack call - dsbtrd 
[D,E,P,INFO] = lapack_dsbtrd(UPLO,A);

%disp(INFO)


% assign output
varargout(1) = {D};
varargout(2) = {E};
varargout(3) = {P};

return

%%
function [D,E,P,INFO] = lapack_dsbtrd(UPLO,A) 
% lapack call 
[esseplus1,n]=size(A);
s=esseplus1-1;
% X is passed to dsbtrd to compute only the first and last s components of
% the reflectors
X=zeros(n,n);
X(1:s,1:s)=speye(s);
X(end-s+1:end,end-s+1:end)=speye(s);
C = lapack('dsbtrd','U',UPLO,n,s, A,esseplus1,zeros(n,1),zeros(n-1,1),X,n,zeros(n,1),0);
[D,E,P,INFO] = C{[7,8,9,12]};
return

