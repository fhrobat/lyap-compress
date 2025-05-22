function xi = poles_Zolotarev(eigmin, eigmax, tol)
%
% Function to compute Zolotarev poles
%


 if eigmin/eigmax > 1e-7
        xi = poles_ellipke(eigmin, eigmax, tol);
 else
        xi = poles_asymptotic(eigmin, eigmax, tol); % may be expensive with respect to poles_ellipke
 end