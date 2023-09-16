function save_experiment_results(exp_data_dir, results, cora_iterates_info)
    assert (nargin >=2, "Must have at least 2 arguments");

    if ~exist(exp_data_dir, 'dir')
        mkdir(exp_data_dir);
    end

    % results is a struct - save it to a mat file
    result_path = exp_data_dir + "/factor_graph_results.mat";
    save(result_path, 'results');
    fprintf('Results saved to %s\n', result_path);

    if (nargin == 3)
        cora_iterates_info_path = exp_data_dir + "/factor_graph_cora_iterates_info.mat";
        save(cora_iterates_info_path, 'cora_iterates_info');
        fprintf('CORA iterates info saved to %s\n', cora_iterates_info_path);
    end
end

