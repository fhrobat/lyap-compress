function xi = poles_asymptotic(a,b, tol)
%
% Function to compute Zolotarev poles using asymptotic expansions 
% 
% Inputs:
%   a,b         the Zolotarev function associated to the poles will be defined over [-b,-a], [a,b] and the poles will belong to [-b,-a]
%   tol         the Zolotarev poles will be computed such that the Zolotarev number is smaller than tol
%
% Outputs:
% xi            vector containing the Zolotarev poles (included in [-b,-a]) 

kp = a/b;
%m_zol = ceil(1/pi^2 * log(4/tol) * log(4/kp));
mu = pi/2 * log(4/kp) / ellipke(kp^2);
m_zol = ceil(1/pi^2 * log(4/tol) * mu);

syms x
K = vpa(subs(log(4/x), kp), 32);
xi = zeros(1,m_zol);

for ii = 1:m_zol
   t = vpa(subs((2*ii-1)*x/(2*m_zol), K) , 32);
   prectanh = vpa(subs(tanh(x), t), 32);
   preccosh = vpa(subs(cosh(x), t), 32);
   precsinh = vpa(subs(sinh(x), 2*t), 32);
   expo = vpa(subs(exp(2*x)/128, t), 32);
   sn = vpa(subs(prectanh + (x^2)/2 * 1/(preccosh^2) * (1/8 * precsinh - 1/4 *t) + (x^2)/8 - (x^4)/4 * expo,kp),32);
   dn = vpa(subs(sqrt(1-x), (1-kp^2/2) * sn^2), 32);
   xi(ii) = vpa(subs(-x*dn, b), 32);
   xi(ii) = double(real(xi(ii)));
   if xi(ii) > -a
       random = rand(1);
       xi(ii) = xi(ii-1) * random - a * (1 - random); % cannot be computed exactly
   end

end

end


