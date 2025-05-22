function param = bcg(A,B,param)
% param = BCG(A,B,param) computes a block conjugate gradients approximation
% to A X = B with parameters specified in the param struct.  Note that the
% classical inner product is hard-corded here; as of Nov 2019, other inner
% products are not yet available.
%
% Possible configurations for param are described in the help line of
% BLOCK_SOLVE. See also BGMRES.
%
% param.bcg_variant can be set to the following possibilities, also known
% as "retooled" bcg variants from
% [Dub2001] A. Dubrulle: Retooling the method of block conjugate gradients,
% Electronic transactions on Numerical Analysis, Vol. 12, pp. 216-233,
% 2001.
%
% The standard algorithm is taken from
% [OLeary1980] D. O'Leary: The block conjugate gradient algorithm and
% related methods, Linear Algebra and its Applications, Vol. 29, pp.
% 293-322, 1980.
%
%   'BCG' (default)         original version proposed in OLeary1980 without 
%                           restarts or stability checks
%   'BCGA'                  first retooled version proposed in Dub2001; a
%                           rearrangement of BCG, but otherwise
%                           mathematically equivalent.
%   'BCGAdQ'                another retooled version proposed in Dub2001,
%                           also referred to therein as BCGAdF with QR; a
%                           slight improvement of BCG by orthogonalization
%                           of the descent matrix
%   'BCGrQ'                 another retooled version proposed in Dub2001;
%                           avoids the need for deflation altogether and
%                           has the best performance in terms of flop
%                           counts and numerical experiments. Based on
%                           orthogonalization of the residual matrix.
%
% Possible configurations for param are described in the help line of
% BLOCKSOLVE.
%
% Further notes:
% We in fact follow the notation of the Ph.D. thesis
% [Birk2015] S. Birk: Deflated shifted block Krylov subspace methods for
% Hermitian positive definite matrices, Fakultät für Mathematik und
% Naturwissenschaften, Bergische Universität Wuppertal, 2015.
% 
% OLeary1980 (Algorithm 6.7 in Birk2015) provides a standard algorithm, but
% it can be unstable due to loss of orthogonality. OLeary1980 recommends
% restarting once rank-deficiency is detected.  We instead implement
% reetooled versions from Dub2001 (see also Algorithms 6.8-6.10 of
% Birk2015).  Versions with deflation are also possible (see Chapter 7 of
% Birk2015).

%%
% Global variables and defaults
global Acount matvecs max_iterations
if isempty(Acount)
    Acount = 0;
end
if isempty(matvecs)
    matvecs = 0;
end
if isfield(param,'max_iterations')
    m = param.max_iterations;                                               % for independent runs
else
    m = max_iterations;                                                     % for use within EXTENDED_KRYLOV
end

if ~isfield(param,'bcg_variant')
    param.bcg_variant = 'bcgrq';
end
param.bcg_variant = lower(param.bcg_variant);                               % allows for typos

if param.verbose
    fprintf('BCG::%s\n',param.bcg_variant);
end

% Initialize
Xm = 0*B;
R = B;
res_norm = 1;                                                               % first approximation is 0
if ~isfield(param,'res_scale')
    res_scale = block_norm(B, param.norm, A);                                % value by which to scale residual
else
    res_scale = param.res_scale;
end

it_count = 0;
switch param.bcg_variant

    case {'bcg', 'standard'}
        %%% BCG
        Rold = R;
        P = Rold;
        while res_norm(end) >= param.tol_res && it_count <= m
            it_count = it_count + 1;
            Z = struct_mult(A,P);
            Acount = Acount + 1;
            matvecs = matvecs + size(P,2);
            M = ((P'*Z)\Rold')*Rold;
            Xm = Xm + P*M;
            R = Rold - Z*M;
            N = ((Rold'*Rold)\R')*R;
            P = R + P*N;
            Rold = R;
            res_norm(end+1) = block_norm(R, param.norm)/res_scale;
            if param.verbose
                fprintf('     it %d: res = %g\n',...
                    it_count, res_norm(end));
            end
        end
    case {'bcga'}
        %%% BCGA
        P = R;
        while res_norm(end) >= param.tol_res && it_count <= m
            it_count = it_count + 1;
            Z = struct_mult(A,P);
            Acount = Acount + 1;
            matvecs = matvecs + size(P,2);
            M = ((P'*Z)\P')*R;
            Xm = Xm + P*M;
            R = R - Z*M;
            N = -((P'*Z)\Z')*R;
            P = R + P *N;
            res_norm(end+1) = block_norm(R, param.norm)/res_scale;
            if param.verbose
                fprintf('     it %d: res = %g\n',...
                    it_count, res_norm(end));
            end
        end
    case {'bcgadf','bcgadq'}
        %%% BCGAdF with QR / BCGAdQ
        [F,C] = qr(R,0);
        while res_norm(end) >= param.tol_res && it_count <= m
            it_count = it_count + 1;
            Z = struct_mult(A,F);
            Acount = Acount + 1;
            matvecs = matvecs + size(F,2);
            M = ((F'*Z)\F')*R;
            Xm = Xm + F*M;
            R = R - Z*M;
            N = -((F'*Z)\Z')*R;
            [F,C] = qr(R + F*N,0);
            res_norm(end+1) = block_norm(C, param.norm)/res_scale;
            if param.verbose
                fprintf('     it %d: res = %g\n',...
                    it_count, res_norm(end));
            end
        end
    case {'bcgrq'}
        %%% BCGrQ
        [Q,C] = qr(R,0);
        D = Q;
        while res_norm(end) >= param.tol_res && it_count <= m
            it_count = it_count + 1;
            Z = struct_mult(A,D);
            Acount = Acount + 1;
            matvecs = matvecs + size(D,2);
            Xm = Xm + (D/(D'*Z))*C;
            [Q,S] = qr(Q - Z/(D'*Z),0);
            D = Q + D*S';
            C = S*C;
            res_norm(end+1) = block_norm(C, param.norm)/res_scale;
            if param.verbose
                fprintf('     it %d: res = %g\n',...
                    it_count, res_norm(end));
            end
        end
end

param.Acount = Acount;
param.matvecs = matvecs;
param.res_norm = res_norm;
param.Xm = Xm;
end