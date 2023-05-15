function save_experiment_results(results, cora_iterates_info, exp_data_path)

    % results is a struct - save it to a mat file
    result_path = strrep(exp_data_path, '.mat', '_results.mat');
    save(result_path, 'results');
    fprintf('Results saved to %s\n', result_path);

    cora_iterates_info_path = strrep(exp_data_path, '.mat', '_cora_iterates_info.mat');
    save(cora_iterates_info_path, 'cora_iterates_info');
    fprintf('CORA iterates info saved to %s\n', cora_iterates_info_path);

end

