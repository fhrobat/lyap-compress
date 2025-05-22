function [param, modified] = param_init_mat_eq(param)
% [param, modified] = PARAM_INIT_MAT_EQ(param) generates and checks the
% input parameter struct.  Fields are checked alphabetically unless a
% logical reason calls for an exception.
%
% param = PARAM_INIT_MAT_EQ returns a default setting
%
% [param, modified] = PARAM_INIT_MAT_EQ(param) returns a corrected parameter
% setting

%%

modified = 0;

if ~nargin
    param = struct;
    modified = 1;
end

if ~isfield(param,'verbose')
    param.verbose = 0;
    disp('Warning: .verbose not specified, set to 1.');
    modified = 1;
end
    
if ~isfield(param,'max_restarts')
    param.max_restarts = 100;
    if param.verbose
        disp('Warning: .max_restarts not specified, set to 100.');
    end
    modified = 1;
end

if ~isfield(param,'memory_max')
    param.memory_max = 100;
    if param.verbose
        disp('Warning: .memory_max not specified, set to 100.');
    end
    modified = 1;
end

if ~isfield(param,'norm')
    param.norm = 'fro';
    if param.verbose
        disp('Warning: .norm not specified, set to ''fro''.');
    end
    modified = 1;
end

if ~isfield(param,'tol_res')
    param.tol_res = 1e-8;
    if param.verbose
        disp('Warning: .tol_res not specified, set to 1e-8.');
    end
    modified = 1;
end

if ~isfield(param,'tol_comp')
    param.tol_comp = param.tol_res*1e-4;
    if param.verbose
        disp('Warning: .tol_comp not specified, set to 1e-4*tol_res.');
    end
    modified = 1;
end

end