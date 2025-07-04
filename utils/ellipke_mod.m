function [k,e] = ellipke_mod(m,tol)
%
% Modified ellipke to accept as input sqrt(1-m) rather then m

if nargin<1
  error(message('MATLAB:ellipke:NotEnoughInputs')); 
end

classin = superiorfloat(m);

if nargin<2, tol = eps(classin); end
if ~isreal(m) || ~isreal(tol)
    error(message('MATLAB:ellipke:ComplexInputs'))
end
if isempty(m), k = zeros(size(m),classin); e = k; return, end
if any(m(:) < 0) || any(m(:) > 1)
  error(message('MATLAB:ellipke:MOutOfRange'));
end
if ~isscalar(tol) || tol < 0 || ~isfinite(tol)
  error(message('MATLAB:ellipke:NegativeTolerance'));
end

a0 = 1;
b0 = m;
c0 = NaN;
s0 = 1-m^2;
i1 = 0; mm = Inf;
while mm > tol
    a1 = (a0+b0)/2;
    b1 = sqrt(a0.*b0);
    c1 = (a0-b0)/2;
    i1 = i1 + 1;
    w1 = 2^i1*c1.^2;
    mm = max(w1(:));
    
    % test for stagnation (may happen for TOL < machine precision)
    if isequal(c0, c1)
        error(message('MATLAB:ellipke:FailedConvergence'));
    end
    
    s0 = s0 + w1;  
    a0 = a1;  b0 = b1;  c0 = c1;
end
k = pi./(2*a1);
e = k.*(1-s0/2);
im = find(m==1);
if ~isempty(im)
    e(im) = ones(length(im),1);
    k(im) = inf;
end
