function build_rail_matrices(src_dir, dst_dir, sizes)
%
% Builds the files dst_dir/rail_N.mat (variables A, E, C) used by
% test_fenics_rail.m from the FEniCS Rail data available at
%       https://zenodo.org/records/5113560
% (download and unzip fenics_rail_N.zip into src_dir).
%
% As in getrail.m (provided with the data):
%   A = -(S + M_GAMMA)   (negative definite),   E = M,
%   C = output matrix (6 x N); its first row is used as right-hand side.
%
% Inputs:
%   src_dir     folder containing the files ODE_N_S.mtx, ODE_N_M_GAMMA.mtx, ODE_N_M.mtx
%   dst_dir     output folder (e.g. the folder 'matrices' of this repository)
%   sizes       sizes N to convert (default [79841, 317377, 1265537])

if nargin < 3
    sizes = [79841, 317377, 1265537];
end
if ~isfolder(dst_dir), mkdir(dst_dir); end

for N = sizes
    base = fullfile(src_dir, sprintf('ODE_%d', N));
    S  = read_mtx([base '_S.mtx']);
    Mg = read_mtx([base '_M_GAMMA.mtx']);
    E  = read_mtx([base '_M.mtx']);
    A  = -(S + Mg);

    C = zeros(6, N);
    C(1, [4, 22, 60])      = [-1, -1, 3];
    C(2, [2, 3, 63])       = [-1, -1, 2];
    C(3, [43, 51])         = [-1, 1];
    C(4, [47, 55])         = [-1, 1];
    C(5, [9, 16, 92])      = [-1, -1, 2];
    C(6, [10, 15, 34, 83]) = [-1, -1, -1, 3];

    save(fullfile(dst_dir, sprintf('rail_%d.mat', N)), 'A', 'E', 'C', '-v7.3');
    fprintf('saved rail_%d.mat\n', N);
end
end


function A = read_mtx(file)
% reads a sparse matrix in MatrixMarket coordinate format (real general)
fid = fopen(file, 'r');
if fid < 0, error('Cannot open %s', file); end
header = fgetl(fid);
if ~contains(header, 'coordinate') || ~contains(header, 'general')
    fclose(fid);
    error('Unsupported MatrixMarket format in %s', file);
end
line = fgetl(fid);
while startsWith(line, '%')
    line = fgetl(fid);
end
sz = sscanf(line, '%d');
data = fscanf(fid, '%d %d %f', [3, sz(3)]);
fclose(fid);
A = sparse(data(1,:), data(2,:), data(3,:), sz(1), sz(2));
end
