function save_experiment_results(results, manopt_info, exp_data_path)

    % results is a struct - save it to a mat file
    result_path = strrep(exp_data_path, '.mat', '_results.mat');
    save(result_path, 'results');
    fprintf('Results saved to %s\n', result_path);

    solver_info_path = strrep(exp_data_path, ".mat", "_solver_info.mat");
    save(solver_info_path, 'manopt_info');
    fprintf('Manopt info saved to %s\n', solver_info_path);

end

