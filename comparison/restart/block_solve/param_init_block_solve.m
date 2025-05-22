function [param, modified] = param_init_block_solve(param)
% [param,modified] = PARAM_INIT_BLOCK_SOLVE(param) generates and checks
% the input parameter struct.
%
% param = PARAM_INIT_BLOCK_SOLVE returns a default setting
%
% param = PARAM_INIT_BLOCK_SOLVE(param) returns a corrected parameter
% setting

%%

modified = 0;

if ~nargin
    param = struct;
    modified = 1;
end

if ~isfield(param,'verbose')
    param.verbose = 1;
    disp('Warning: .verbose not specified, set to 1.');
    modified = 1;
end

if ~isfield(param,'max_iterations')
    param.max_iterations = 100;
    if param.verbose
        disp('Warning: .max_iterations not specified, set to 100.');
    end
end

if ~isfield(param,'max_restarts')
    param.max_restarts = 100;
    if param.verbose
        disp('Warning: .max_restarts not specified, set to 100.');
    end
end

if ~isfield(param,'norm')
    param.norm = 'fro';
    if param.verbose
        disp('Warning: .norm not specified, set to ''fro''.');
    end
end

if ~isfield(param,'solver')
    param.solver = 'bgmres';
    disp('Warning: .solver not specified, set to ''bgmres''.')
end

if ~isfield(param,'tol_res')
    param.tol_res = 1e-8;
    if param.verbose
        disp('Warning: .tol_res not specified, set to 1e-8.');
    end
    modified = 1;
end

if ~isfield(param,'tol_zero')
    param.tol_zero = 1e-13;
    if param.verbose
        disp('Warning: .tol_zero not specified, set to 1e-13.');
    end
    modified = 1;
end
end