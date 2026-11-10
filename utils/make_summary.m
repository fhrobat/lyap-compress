function make_summary(out_dir)
%
% Collects the results saved by the test scripts in out_dir/raw (one file per
% experiment and size) into one csv file per experiment, out_dir/<test>.csv,
% and writes out_dir/summary.md with:
%   - the comparison of compress, two-pass and restart: number of
%     matrix-vector products, time, time ratio compress / two-pass, residual;
%   - the results of all the methods;
%   - the sensitivity experiment;
%   - the exit status of every MATLAB run (from out_dir/jobs.csv, if present).

raw = fullfile(out_dir, 'raw');
tests  = {'test_2DLap', 'test_3DLap', 'test_fenics_rail', 'test_heat'};
titles = {'2D Laplacian', '3D Laplacian', ...
          'FEniCS Rail', 'Heat equation'};

fid = fopen(fullfile(out_dir, 'summary.md'), 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '# Results\n\n');
fprintf(fid, 'Folder: `%s`. Machine and MATLAB version: see `machine.txt`.\n\n', out_dir);
fprintf(fid, 'matvecs of two-pass = matrix-vector products of both passes; ');
fprintf(fid, 'residual = scaled residual norm computed a posteriori; ratio = time compress / time two-pass.\n\n');

for t = 1:numel(tests)
    fprintf(fid, '## %s\n\n', titles{t});
    T = collect(raw, tests{t});
    if isempty(T)
        fprintf(fid, 'No results.\n\n');
        continue
    end
    writetable(T, fullfile(out_dir, [tests{t} '.csv']));

    % compress, two-pass and restart
    fprintf(fid, '| N | k | matvecs compress | matvecs two-pass | matvecs restart | time compress | time two-pass | time restart | ratio | residual compress | residual two-pass |\n');
    fprintf(fid, '|---|---|---|---|---|---|---|---|---|---|---|\n');
    for N = unique(T.N, 'stable')'
        tc = val(T, N, "compress", 'time');
        tt = val(T, N, "two-pass", 'time');
        fprintf(fid, '| %d | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n', N, ...
            num(val(T, N, "compress", 'num_poles'), '%d'), ...
            num(val(T, N, "compress", 'matvecs'), '%d'), num(val(T, N, "two-pass", 'matvecs'), '%d'), ...
            num(val(T, N, "restart", 'matvecs'), '%d'), ...
            num(tc, '%.1f'), num(tt, '%.1f'), num(val(T, N, "restart", 'time'), '%.1f'), ...
            num(tc/tt, '%.2f'), ...
            num(val(T, N, "compress", 'true_res'), '%.1e'), num(val(T, N, "two-pass", 'true_res'), '%.1e'));
    end
    fprintf(fid, '\n');

    % all the methods
    fprintf(fid, '<details><summary>All methods</summary>\n\n');
    fprintf(fid, '| N | method | time | matvecs | residual | poles |\n|---|---|---|---|---|---|\n');
    for i = 1:height(T)
        fprintf(fid, '| %d | %s | %s | %s | %s | %s |\n', T.N(i), T.method(i), num(T.time(i), '%.1f'), ...
            num(T.matvecs(i), '%d'), num(T.true_res(i), '%.1e'), num(T.num_poles(i), '%d'));
    end
    fprintf(fid, '\n</details>\n\n');
end

% sensitivity experiment
fprintf(fid, '## Sensitivity to the estimate of the smallest eigenvalue\n\n');
f = fullfile(raw, 'test_2DLap_sensitivity.csv');
if isfile(f)
    S = readtable(f, 'TextType', 'string');
    copyfile(f, fullfile(out_dir, 'test_2DLap_sensitivity.csv'));
    names = S.Properties.VariableNames;
    fprintf(fid, '| %s |\n', strjoin(names, ' | '));
    fprintf(fid, '|%s\n', repmat('---|', 1, numel(names)));
    for i = 1:height(S)
        cells = cell(1, numel(names));
        for jj = 1:numel(names)
            x = S.(names{jj})(i);
            if isstring(x), cells{jj} = char(x);
            elseif x == round(x), cells{jj} = sprintf('%d', x);
            else, cells{jj} = sprintf('%.3g', x);
            end
        end
        fprintf(fid, '| %s |\n', strjoin(cells, ' | '));
    end
    fprintf(fid, '\n');
else
    fprintf(fid, 'No results.\n\n');
end

% exit status of the MATLAB runs
f = fullfile(out_dir, 'jobs.csv');
if isfile(f)
    o = detectImportOptions(f, 'Delimiter', ',');
    o = setvartype(o, {'script', 'size'}, 'string');
    J = readtable(f, o);
    fprintf(fid, '## MATLAB runs\n\n');
    fprintf(fid, 'Exit status 0 = completed; 137 = killed (out of memory, limit in `machine.txt`); other = error, see `logs/`.\n\n');
    fprintf(fid, '| script | size | exit status | duration (s) |\n|---|---|---|---|\n');
    for i = 1:height(J)
        fprintf(fid, '| %s | %s | %d | %d |\n', J.script(i), J.size(i), J.status(i), J.seconds(i));
    end
    fprintf(fid, '\n');
end
end


function T = collect(raw, test)
% concatenates raw/<test>_size<i>.csv in the order of the sizes
files = dir(fullfile(raw, [test '_size*.csv']));
T = table();
if isempty(files), return; end
idx = cellfun(@(s) sscanf(s, [test '_size%d.csv']), {files.name});
[~, p] = sort(idx);
for f = files(p)'
    T = [T; readtable(fullfile(f.folder, f.name), 'TextType', 'string')]; %#ok<AGROW>
end
end

function v = val(T, N, method, field)
% value of field for the given size and method (NaN if missing)
i = find(T.N == N & T.method == method, 1);
if isempty(i), v = NaN; else, v = T.(field)(i); end
end

function s = num(x, fmt)
% formats a number, '-' if missing
if isempty(x) || isnan(x), s = '-'; else, s = sprintf(fmt, x); end
end
