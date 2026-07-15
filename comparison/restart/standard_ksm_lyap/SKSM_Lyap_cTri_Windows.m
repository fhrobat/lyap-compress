function [Z,time_tot,time_res]=SKSM_Lyap_cTri_Windows(A,C,m,tol)
%[Z,time_tot,time_res]=SKSM_Lyap_cTri_Windows(A,C,m,tol)
% Solution of the Lyapunov equation 
%
% A*X+X*A+C*C^T=0,    A sym neg.def. and symmetric, C tall
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
% email: davide.palitta3@unibo.it, valeria.simoncini@unibo.it
%
% REQUIREMENT: the subroutines eig_dstevr and tridiag_dsbtrd which make use of LAPACK and C-BLAS subroutines. 
% 
%
% INPUT
% A: sym negative def. matrix  
% C: rhs matrix 
% m: max number of iterations
% tol: final residual accuracy (relative residual F-norm)
%
% OUTPUT
% Z: tall matrix such that Z*Z^T is the approx solution (Z is truncated to sing.val >1e-12)
% time_tot: total cputime
% time_res: total time requested to compute the residual norm
%
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS 
% FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
% COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER 
% IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 

% Check the input
if norm(A-A','fro')>1e-12 
    error('The coefficient matrix must be symmetric')
end

tic;
[n,s]=size(C); 
nrmb=norm(C,'fro')^2;

% Orthonormalize C
[U(1:n,1:s),beta]=qr(C,0);

beta2=beta*beta';

% Preallocation
TA=zeros((m+1)*s,m*s);
TTA=zeros(s+1,m*s);

time_res=zeros(m,1);
time_subtot=0;
        
t1=toc;
for j=1:m
    tic
    jms=(j-1)*s+1;j1s=(j+1)*s;js=j*s;js1=js+1;
   
    %sequence in A
    %new basis block (double modified gram-schmidt)
    Up = A*U(1:n,end-s+1:end);
  
    for l=1:2
        K=min([j, 3]);
        for kk=1:K
              k1=(kk-1)*s+1; k2=kk*s;
              coeff= U(1:n,k1:k2)'*Up;
              TA((j-K+kk-1)*s+1:(j-K+kk)*s,jms:js) = TA((j-K+kk-1)*s+1:(j-K+kk)*s,jms:js)+ coeff; 
              Up = Up - U(1:n,k1:k2)*coeff; 
        end
    end
    
  
    if (j<2)
        [U(1:n,s+1:2*s),TA(js1:j1s,jms:js)]=qr(Up,0);
        
        % prepare LLA for dstrd
        for i=1:s
            TTA(end-i+1,i:s)=diag(TA(jms:js,jms:js),i-1);
        end
        LLA=TTA(:,1:j*s);
        
    else
        % move the blocks
        if j>2
            U(1:n,1:2*s)=U(1:n,s+1:3*s);
        end
        [U(1:n,2*s+1:3*s),TA(js1:j1s,jms:js)]=qr(Up,0);
    
        % prepare LLA for dstrd
        HHA=zeros(s+1,s);
        for i=1:s+1
            HHA(end-i+1,1:s)=diag(TA(js-s+2-i:js-i+1,js-s+1:js));
        end
        TTA(:,(j-1)*s+1:j*s)=HHA;
        LLA=TTA(:,1:j*s);
    
    end
    
    t2=toc;
        
    tic
    % take the ((j+1)s,js) block of TA 
    ll = TA(js1:j1s,jms:js)';             

    % tridiagonalize LLA by dsbtrd
    [D,E,P]=tridiag_dsbtrd(LLA);
    
    % compute the eigendecompositions of the obtained tridiagonal matrices by dstevr
    [Q,lambda]=eig_dstevr(D,E);
    
    
    % Compute the residual norm
    S=((P(1:s,:)*Q)'*beta2*(P(1:s,:)*Q))./(lambda*ones(1,js)+ones(js,1)*lambda');
    W=(P(end-s+1:end,:)*Q)'*ll;
    
    res=sqrt(2)*norm(S*W,'fro');

    
    time_res(j)=toc;

    fprintf('It: %3d, Relative Res: %10.5e\n', j, res/nrmb)
    time_subtot=time_subtot+t2;
    
    if (res/nrmb<tol) 
         break
    end
    

end
tic
it=j;


% Compute the full eigendecomposition of T and the matrix Y
TA(1:js,1:js)=(TA(1:js,1:js)+TA(1:js,1:js)')/2;
[QA,lambda]=eig(full(TA(1:js,1:js)));
lambdaA=diag(lambda);
Q=QA(1:s,:);

% Compute the solution of the last projected equation
Y=-(Q'*beta2*Q)./(lambdaA*ones(1,js)+ones(js,1)*lambdaA');

% Truncation strategy
% Truncation strategy
[uY,sY]=eig(Y); [sY,id]=sort(abs(diag(sY)));
sY=flipud(sY);
uY=uY(:,id(end:-1:1));
is=sum(abs(sY)>1e-12);
Y0 = QA*(uY(:,1:is)*diag(sqrt(sY(1:is))));
t4=toc;

% Two-pass

U=[];
tic
[U(1:n,1:s),~]=qr(C,0);
Z=U(1:n,1:s)*Y0(1:s,:);

%separe j=1,2 from the other j's to avoid "if" within the "for" loop
j=1;
j1s=2*s;
js1=s+1;
% sequence in A
UpA = A*U(1:n,end-s+1:end);
U(1:n,s+1:2*s)=(UpA-U(1:n,end-s+1:end)*TA((j-1)*s+1:j*s,(j-1)*s+1:j*s))/TA(j*s+1:(j+1)*s,(j-1)*s+1:j*s);
Z=Z+U(:,end-s+1:end)*Y0(js1:j1s,:);

j=2;
j1s=3*s;
js1=2*s+1;
%sequence in A
UpA = A*U(1:n,end-s+1:end);
U(1:n,2*s+1:3*s)=(UpA-U(1:n,1:2*s)*TA((j-2)*s+1:j*s,(j-1)*s+1:j*s))/TA(j*s+1:(j+1)*s,(j-1)*s+1:j*s);
Z=Z+U(:,end-s+1:end)*Y0(js1:j1s,:);

for j=3:it-1
    j1s=(j+1)*s;
    js1=j*s+1;
    
    % sequence in A
    UpA = A*U(1:n,end-s+1:end);
    U(1:n,1:2*s)=U(1:n,s+1:3*s);
    U(1:n,2*s+1:3*s)=(UpA-U(1:n,1:2*s)*TA((j-2)*s+1:j*s,(j-1)*s+1:j*s))/TA(j*s+1:(j+1)*s,(j-1)*s+1:j*s);
    Z=Z+U(:,end-s+1:end)*Y0(js1:j1s,:);
end
 

time_twopass=toc;
fprintf('\n time_twopass: %10.5e\n',time_twopass)

time_subtot=time_subtot+t4+t1+time_twopass;
time_res=sum(time_res);
time_tot=time_subtot+time_res;    