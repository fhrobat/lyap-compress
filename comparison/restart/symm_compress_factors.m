function [U, L] = symm_compress_factors(U, L, param)
% [U, L] = SYMM_COMPRESS_FACTORS(U, L, tol) returns a compressed version of
% U and core L, compressed up to tolerance tol.

%%
[Qu, Ru] = qr(full(U), 0);
[V,  D]  = eig(full(Ru * L * Ru'));
d=diag(D);
[~,index]=sort(abs(d),'descend');
d=d(index);
V=V(:,index);
tol=param.tol_comp;
switch param.norm
        case 'fro'
            dd=flipud(d);
            ind=find(sqrt(cumsum(dd.^2)) > tol,1,'first');
            
            if ~isempty(ind)
                ind=length(d)-ind+1;
            else
                ind=length(d);
            end
           
        case {2, '2'}
            ind = find(abs(d) > tol,1,'last');
            if isempty(ind)
               ind=1;
            end
end 
U  = Qu * V(:, 1:ind);
L = diag(d(1:ind));
L = (L + L')/2;
end