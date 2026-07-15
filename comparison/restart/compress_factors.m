function [U, V] = compress_factors(U, V, param)
% [U, V] = COMPRESS_FACTORS(U, V, tol) returns a compressed version of
% U and core L, compressed up to tolerance tol.

%%
[Qu, Ru] = qr(full(U),0);
[Qv, Rv] = qr(full(V),0);
[W, S, Z] = svd(full(Ru * Rv'));
tol=param.tol_comp;
switch param.norm
        case 'fro'
            dd=flipud(diag(S));
            ind=find(sqrt(cumsum(dd.^2)) > tol,1,'first');
            
            if ~isempty(ind)
                ind=size(S,1)-ind+1;
            else
                ind=size(S,1);
            end
           
        case {2, '2'}
            ind = find(diag(S) > tol,1,'last');
            if isempty(ind)
               ind=1;
            end
end 
r = sum(diag(S > tol));
U = Qu * W(:, 1:ind) * sqrt(S(1:ind, 1:ind));
V = Qv * Z(:, 1:ind) * sqrt(S(1:ind, 1:ind));
end