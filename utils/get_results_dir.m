function results_dir = get_results_dir(root)
%
% Folder where the test scripts save their results: the environment variable
% LYAP_RESULTS_DIR if it is set (run_all.sh sets it to results/<machine>_<date>/raw),
% otherwise root/results.

results_dir = getenv('LYAP_RESULTS_DIR');
if isempty(results_dir)
    results_dir = fullfile(root, 'results');
end
if ~isfolder(results_dir)
    mkdir(results_dir);
end
end
