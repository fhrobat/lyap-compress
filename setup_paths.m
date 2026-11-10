function root = setup_paths()
%
% Adds to the MATLAB path the folders needed by lyap_compress and by the
% test scripts.
%
% The compress-and-restart code (https://gitlab.com/katlund/compress-and-restart-KSM),
% needed for the methods restart and extended, is expected in
%   external/compress-and-restart-KSM
% (see README.md).

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'lyap_compress'));
addpath(genpath(fullfile(root, 'utils')));
addpath(genpath(fullfile(root, 'comparison')));

ext = fullfile(root, 'external', 'compress-and-restart-KSM', 'main');
if isfolder(ext)
    addpath(genpath(ext));
end
end
