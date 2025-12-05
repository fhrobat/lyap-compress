function [U, nmatvecs] = rational_krylov_pcg(A, W, poles, options)
	% constructs rational Krylov subspace Q(A, W, poles) solving linear systems employing pcg with ichol preconditioning
    % tol indicates the tolerance fed into pcg
	% by default the starting vector W is not included (i.e., no starting infinite pole is added) 
	
    
    maxit = 1000;
    nmatvecs = 0;


	if nargin < 4
		options = [];
    end

    if isfield(options, 'tol_solver') == 0
        options.tol_solver = 1e-6;
    end


	if isfield (options, 'isreal') == 0
		options.isreal = 0;
	end

	if abs(poles(1)) == inf
		W = W/ norm(W);
		j = 2;
		U = W;
	elseif options.isreal && poles(2) == conj(poles(1))
		W = W/ norm(W);
		%W = (A - poles(1)*speye(size(A))) \  W;
        %%
        if strcmp(options.precond, 'ichol')
            L = ichol(A - poles(1)*speye(size(A)));
            [W,~,~,it] = pcg(A - poles(1)*speye(size(A)),W,options.tol_solver,maxit,L,L');
        else
            [W,~,~,it] = pcg(A - poles(1)*speye(size(A)),W,options.tol_solver,maxit);
        end
        nmatvecs = nmatvecs + it;
        %%
		W1=real(W);
		W=imag(W);
		j = 3;
		[U,~] = qr([W1, W], 0);
		W=U(:,end);
	else
		W = W/ norm(W);
		%W = (A - poles(1)*speye(size(A))) \  W;
        %%
        if strcmp(options.precond, 'ichol')
            L = ichol(A - poles(1)*speye(size(A)));
            [W,~,~,it] = pcg(A - poles(1)*speye(size(A)),W,options.tol_solver,maxit,L,L');
        else
            [W,~,~,it] = pcg(A - poles(1)*speye(size(A)),W,options.tol_solver,maxit);
        end
        nmatvecs = nmatvecs + it;
        %%
		j = 2;
		W = W/ norm(W);
		U = W;
	end

	while j <= length(poles)
		if abs(poles(j)) == inf
			W = A * W;
			j = j+1;
		elseif options.isreal && j+1<=length(poles) && poles(j+1) == conj(poles(j))
			%W = (A - poles(j)*speye(size(A))) \  W;
            %%
            if strcmp(options.precond, 'ichol')
                L = ichol(A - poles(j)*speye(size(A)));
                [W,~,~,it] = pcg(A - poles(j)*speye(size(A)),W,options.tol_solver,maxit,L,L');
            else
                [W,~,~,it] = pcg(A - poles(j)*speye(size(A)),W,options.tol_solver,maxit);
            end
            nmatvecs = nmatvecs + it;
            %%
			W1 = real(W);
			W = imag(W);
			
			% cgs2
			W1 = W1 - U * (U'*W1);
			W1 = W1 - U * (U'*W1);
			W1 = W1/ norm(W1);
			U = [U, W1];

			% qr alternative
			% [U,~] = qr([U, W1], 0);
			j = j+2;
		else
			%W = (A - poles(j)*speye(size(A))) \  W;
            %%
            if strcmp(options.precond, 'ichol')
                L = ichol(A - poles(j)*speye(size(A)));
                [W,~,~,it] = pcg(A - poles(j)*speye(size(A)),W,options.tol_solver,maxit,L,L');
            else
                [W,~,~,it] = pcg(A - poles(j)*speye(size(A)),W,options.tol_solver,maxit);
            end
            nmatvecs = nmatvecs + it;
            %%
			j = j+1;
		end

		% cgs2
		W = W - U * (U'*W);
		W = W - U * (U'*W);
		W = W/norm(W);
		U = [U, W];

		% qr alternative
		% [U,~] = qr([U, W], 0);

		W = U(:,end);
    end
end