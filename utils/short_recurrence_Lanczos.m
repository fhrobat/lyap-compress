function [v3,alpha2,beta2] = short_recurrence_Lanczos(mult, v1, v2, beta1, beta2, alpha2)
%
% Runs one iteration of Lanczos algorithm employing short recurrence
%
% Inputs:
%	 mult			function handle such that mult(v) = A*v
%	 v1, v2			last two basis vectors
%	 beta1			current last entry of beta vector
%    beta2          next beta coefficient (optional)
%    alpha2         next alpha coefficient (optional)
%
%
% Outputs:
%	v3				next basis vector
%	alpha2, beta2	next entries in alpha and beta vectors

if nargin < 5
    tildey = mult(v2)-beta1*v1;
    alpha2 = (tildey'*v2);
    w = tildey-alpha2*v2;
    beta2 = norm(w,2);
    v3 = w / beta2;
else
    tildey = mult(v2)-beta1*v1;
    w = tildey-alpha2*v2;
    v3 = w / beta2;
end
end