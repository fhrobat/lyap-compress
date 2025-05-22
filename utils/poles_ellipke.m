function xi = poles_ellipke(a,b,tol)
%
% Function to compute Zolotarev poles using MATLAB functions ellipke and ellipj
% 
% Inputs:
%   a,b         the Zolotarev function associated to the poles will be
%   defined over [-b,-a], [a,b] and the poles will belong to [-b,-a]
%   tol         the Zolotarev poles will be computed such that the Zolotarev number is smaller than tol
%
% Outputs:
% xi            vector containing the Zolotarev poles (included in [-b,-a])


kp = a/b;
%m_zol = ceil(1/pi^2 * log(4/tol) * log(4/kp));
K = ellipke_mod(kp);
mu = pi/2 * K/ellipke(kp^2);
m_zol = ceil(1/pi^2 * log(4/tol) * mu);
xi = zeros(1,m_zol);

for ii=1:m_zol
    [~,~,dn]=ellipj_mod((2*ii-1)*K/(2*m_zol),kp);
    xi(ii) = -b*dn;
end


end