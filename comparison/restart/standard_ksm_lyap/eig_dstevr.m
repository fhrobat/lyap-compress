function [varargout] = eig_dstevr(varargin)
%  DSTEVR computes selected eigenvalues and, optionally, eigenvectors
%  of a real symmetric tridiagonal matrix T.  Eigenvalues and
%  eigenvectors can be selected by specifying either a range of values
%  or a range of indices for the desired eigenvalues.
% 
%  Whenever possible, DSTEVR calls DSTEMR to compute the
%  eigenspectrum using Relatively Robust Representations.  DSTEMR
%  computes eigenvalues by the dqds algorithm, while orthogonal
%  eigenvectors are computed from various "good" L D L^T representations
%  (also known as Relatively Robust Representations). Gram-Schmidt
%  orthogonalization is avoided as far as possible. More specifically,
%  the various steps of the algorithm are as follows. For the i-th
%  unreduced block of T,
%     (a) Compute T - sigma_i = L_i D_i L_i^T, such that L_i D_i L_i^T
%          is a relatively robust representation,
%     (b) Compute the eigenvalues, lambda_j, of L_i D_i L_i^T to high
%         relative accuracy by the dqds algorithm,
%     (c) If there is a cluster of close eigenvalues, "choose" sigma_i
%         close to the cluster, and go to step (a),
%     (d) Given the approximate eigenvalue lambda_j of L_i D_i L_i^T,
%         compute the corresponding eigenvector by forming a
%         rank-revealing twisted factorization.
%  The desired accuracy of the output can be specified by the input
%  parameter ABSTOL.
% 
%  For more details, see "A new O(n^2) algorithm for the symmetric
%  tridiagonal eigenvalue/eigenvector problem", by Inderjit Dhillon,
%  Computer Science Division Technical Report No. UCB//CSD-97-971,
%  UC Berkeley, May 1997.
% 
% 
%  Note 1 : DSTEVR calls DSTEMR when the full spectrum is requested
%  on machines which conform to the ieee-754 floating point standard.
%  DSTEVR calls DSTEBZ and DSTEIN on non-ieee machines and
%  when partial spectrum requests are made.
% 
%  Normal execution of DSTEMR may create NaNs and infinities and
%  hence may abort due to a floating point exception in environments
%  which do not handle NaNs and infinities in the ieee standard default
%  manner.

if nargin<1,
	varargout = {};
	help eig_dtevr;
	return;
end
D=varargin{1};
E=varargin{2};
% lapack call - dstevr 
[Q,lambda] = lapack_dstevr(D,E);

%disp(INFO)


% assign output
varargout(1) = {Q};
varargout(2) = {lambda};

return

%%
function [Q,lambda] = lapack_dstevr(D,E) 
% lapack call 
n=length(D);
C = lapack('dstevr','V','A',n,D,E,1,n,1,n,0,0,zeros(n,1),zeros(n,n),n,zeros(2*n,1),zeros(n,1),-1,zeros(n,1),-1,0);
WORK = C{16}; % find optimal size
LWORK = WORK(1);
IWORK=C{18};
LIWORK=IWORK(1);
C = lapack('dstevr','V','A',n,D,E,1,n,1,n,0,0,zeros(n,1),zeros(n,n),n,zeros(2*n,1),zeros(LWORK,1),LWORK,zeros(LIWORK,1),LIWORK,0);
[Q,lambda] = C{[13,12]};
return

% 