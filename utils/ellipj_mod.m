function [sn,cn,dn] = ellipj_mod(u,m,tol)
%
% Modified ellipj to accept as input sqrt(1-m) rather then m


if nargin<2
  error(message('MATLAB:ellipj:NotEnoughInputs')); 
end

classin = superiorfloat(u,m);

if nargin<3, tol = eps(classin); end


if ~isreal(u) || ~isreal(m) || ~isreal(tol)
    error(message('MATLAB:ellipj:ComplexInputs'))
end

if isscalar(m), m = m(ones(size(u))); end
if isscalar(u), u = u(ones(size(m))); end
if ~isequal(size(m),size(u)) 
  error(message('MATLAB:ellipj:InputSizeMismatch')); 
end

mmax = numel(u);

cn = zeros(size(u),classin);
sn = cn;
dn = sn;
m = m(:).';    % make a row vector
u = u(:).';

% check for ~all(m >= 0) instead of m < 0 so that m = NaN throws error
if ~all(m >= 0) || any(m > 1)
  error(message('MATLAB:ellipj:MOutOfRange'));
end
if ~isscalar(tol) || tol<0 || ~isfinite(tol)
  error(message('MATLAB:ellipj:NegativeTolerance'));
end

% pre-allocate space and augment if needed
chunk = 10;
a = zeros(chunk,mmax);
c = a;
b = a;
a(1,:) = ones(1,mmax);
if m > 1e-4
    c(1,:) = sqrt(1-m^2);
else
    c(1,:) = 1 - m^2/2;
end
b(1,:) = m;
n = zeros(1,mmax);
i = 1;
while any(abs(c(i,:)) > tol)
    i = i + 1;
    if i > size(a,1)
      a = [a; zeros(chunk,mmax)]; %#ok<AGROW>
      b = [b; zeros(chunk,mmax)]; %#ok<AGROW>
      c = [c; zeros(chunk,mmax)]; %#ok<AGROW>
    end
    a(i,:) = 0.5 * (a(i-1,:) + b(i-1,:));
    b(i,:) = sqrt(a(i-1,:) .* b(i-1,:));
    c(i,:) = 0.5 * (a(i-1,:) - b(i-1,:));
    
    % test for stagnation (may happen for TOL < machine precision)
    if isequal(c(i, :), c(i-1, :))
        error(message('MATLAB:ellipj:FailedConvergence'));
    end
    
    in = find((abs(c(i,:)) <= tol) & (abs(c(i-1,:)) > tol));
    if ~isempty(in)
      [mi,ni] = size(in);
      n(in) = repmat((i-1), mi, ni);
    end
end
phin = zeros(i,mmax,classin);
phin(i,:) = (2 .^ n).*a(i,:).*u;
while i > 1
    i = i - 1;
    in = find(n >= i);
    phin(i,:) = phin(i+1,:);
    if ~isempty(in)
      phin(i,in) = 0.5 * ...
      (asin(c(i+1,in).*sin(rem(phin(i+1,in),2*pi))./a(i+1,in)) + phin(i+1,in));
    end
end
sn(:) = sin(rem(phin(1,:),2*pi));
cn(:) = cos(rem(phin(1,:),2*pi));
dn(:) = sqrt(1 - (1-m^2) .* (sn(:).').^2);

% special case m = 1 
m1 = find(m==1);
sn(m1) = tanh(u(m1));
cn(m1) = sech(u(m1));
dn(m1) = sech(u(m1));
% special case m = 0
dn(m==0) = 1;
end