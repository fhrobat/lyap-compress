function [U, L] = symm_compress_clean(U, L)
% [U, L] = SYMM_COMPRESS_CLEAN(U, L) returns a compressed version of U and
% core L, by discarding the eigenvectors related to the negative
% eigenvalues.

%%
[Qu, Ru] = qr(full(U), 0);
[V,  D]  = eig(full(Ru * L * Ru'));

% Make sure everything is in the right order
[d,index] = sort(diag(D), 'descend');
V = V(:,index);

% Take the index related to the last nonnegative eigenvalues
ind = find(d>=0,1,'last');
U  = Qu * V(:, 1:ind);
L = diag(d(1:ind));

end