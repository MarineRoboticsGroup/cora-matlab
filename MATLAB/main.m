clc; clear all; close all;
base_data_dir = "~/experimental_data/plaza/Plaza2";

manopt_opts.init = "random";
manopt_opts.verbosity = 2;
manopt_opts.debug = 0;
manopt_opts = get_manopt_opts(manopt_opts);
do_not_lift = false; % whether or not to perform lifting (e.g. vanilla solve vs CORA)
animation_show_gt = true;
plot_show_gt = true;
look_for_cached_soln = true;

% set init strategy and check if it is valid
assert(manopt_opts.init == "random" || manopt_opts.init == "odom" || manopt_opts.init == "gt");
experiments = get_experiment_filepaths_in_subdirs(base_data_dir);
num_experiments = length(experiments);

for exp_idx = 1:num_experiments
    fprintf("Starting experiment %d of %d\n", exp_idx, num_experiments);
    data_path = experiments(exp_idx);
    res_path = strrep(data_path, ".mat", "_results.mat");
    cora_iterates_info_path = strrep(data_path, '.mat', '_cora_iterates_info.mat');
    data_dir = fileparts(data_path);
    exp_data = load_ra_slam_problem(data_path);

    if look_for_cached_soln && exist(res_path, 'file') && exist(cora_iterates_info_path, 'file')
        warning("Results already found for experiment %s ... skipping", data_path);

        % write the obtained result to .tum format for comparison
        res = load(res_path);
        X = res.results.X;
        write_result_to_tum(X, exp_data, data_dir)
        animate_lifted_solver_trajectory(data_path, animation_show_gt);
        continue
    end

    [X, X_is_optimal, cora_iterates_info, Manopt_opts] = ra_slam(exp_data, manopt_opts, do_not_lift);
    res = struct();

    % values from the true final solution
    res.X = X; % the final solution (rounded and refined)
    res.X_is_certified = X_is_optimal;

    % save the results and write to .tum format for comparison
    save_experiment_results(res, cora_iterates_info, data_path);
    write_result_to_tum(X, exp_data, data_dir);

    % visualize the solution
    animate_lifted_solver_trajectory(data_path, animation_show_gt);
end
