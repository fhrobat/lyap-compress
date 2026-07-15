function print_flag(flag, fID)
% PRINT_FLAG(flag) translates flag into a convergence statement.

%%
if nargin == 1
    switch flag
        case 0
            fprintf('Flag unchanged.  Check code for bugs.\n');
        case 1
            fprintf('Method converged to desired tolerance.\n');
        case 2
            fprintf('Maximum number of iterations reached before convergence.\n');
    end
elseif nargin == 2
    switch flag
        case 0
            fprintf(fID, 'Flag unchanged.  Check code for bugs.\n');
        case 1
            fprintf(fID, 'Method converged to desired tolerance.\n');
        case 2
            fprintf(fID, 'Maximum number of iterations or restarts reached before convergence.\n');
    end
end
end