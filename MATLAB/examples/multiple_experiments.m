clc; clear all; close all;

% the data directory is two levels up from this file, let's get its path
central_data_dir = fileparts(fileparts(mfilename('fullpath')));

% the experiments we want are nested under <central_data_dir>/data/MATLAB/simulated/100loop_closures
base_data_dir = fullfile(central_data_dir, 'data', 'MATLAB', 'simulated', '100loop_closures');

% we will set solver params outside of the loop
manopt_opts.init = "random";
manopt_opts.verbosity = 2;
manopt_opts.debug = 0;
manopt_opts = get_manopt_opts(manopt_opts);
assert(manopt_opts.init == "random" || manopt_opts.init == "odom" || manopt_opts.init == "gt");

% set visualization params
animation_show_gt = true;
look_for_cached_soln = false; % if true, will not re-solve problems for which results already exist

% get the list of experiments
experiments = get_experiment_filepaths_in_subdirs(base_data_dir);
num_experiments = length(experiments);

for exp_idx = 1:num_experiments
    fprintf("Starting experiment %d of %d\n", exp_idx, num_experiments);
    data_path = experiments(exp_idx);
    res_path = strrep(data_path, ".mat", "_results.mat");
    cora_iterates_info_path = strrep(data_path, '.mat', '_cora_iterates_info.mat');
    data_dir = fileparts(data_path);
    exp_data = load_ra_slam_problem(data_path);

    % check for cached results
    if look_for_cached_soln && exist(res_path, 'file') && exist(cora_iterates_info_path, 'file')
        warning("Results already found for experiment %s ... skipping", data_path);
        res = load(res_path);
        X = res.results.X;
        write_result_to_tum(X, exp_data, data_dir)
        animate_lifted_solver_trajectory(data_path, animation_show_gt);
        continue
    end

    % run the solver
    [X, X_is_optimal, cora_iterates_info, Manopt_opts] = cora(exp_data, manopt_opts);

    % values from the true final solution
    res = struct();
    res.X = X; % the final solution (rounded and refined)
    res.X_is_certified = X_is_optimal;

    % save the results and write to .tum format for comparison
    save_experiment_results(res, cora_iterates_info, data_path);
    write_result_to_tum(X, exp_data, data_dir);

    % visualize the solution
    animate_lifted_solver_trajectory(data_path, animation_show_gt);
end
