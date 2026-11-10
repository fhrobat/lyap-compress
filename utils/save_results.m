function save_results(T, results_dir, name)
%
% Saves the table T in results_dir/name.csv and results_dir/name.mat.
% The test scripts call it after every method, so that the results obtained
% so far are not lost if MATLAB is stopped (e.g. because it runs out of memory).

writetable(T, fullfile(results_dir, [name '.csv']));
save(fullfile(results_dir, [name '.mat']), 'T');
end
