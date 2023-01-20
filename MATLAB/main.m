clc; clear all; close all;
base_data_dir = "~/data/manhattan/cert/no_loop_closures";
% base_data_dir = "~/data/manhattan/cert/100loop_closures";
% base_data_dir = "~/data/hat_data/16OCT2022";

manopt_opts.init = "odom";
manopt_opts.verbosity = 2;
manopt_opts = get_manopt_opts(manopt_opts);
do_not_lift = false; % whether or not to perform lifting (e.g. vanilla solve vs CORA)


% set init strategy and check if it is valid
assert(manopt_opts.init == "random" || manopt_opts.init == "odom" || manopt_opts.init == "gt");
experiments = get_experiment_filepaths_in_subdirs(base_data_dir);
num_experiments = length(experiments);


for exp_idx = 1:num_experiments
    fprintf("Starting experiment %d of %d\n", exp_idx, num_experiments);
    data_path = experiments(exp_idx);
    res_path = strrep(data_path, ".mat", "_results.mat");
    solver_info_path = strrep(data_path, ".mat", "_solver_info.mat");
    if exist(res_path, 'file') && exist(solver_info_path, 'file')
        warning("Results already found for experiment %s ... skipping", data_path);
        animate_lifted_solver_trajectory(data_path);
        continue
    end

    exp_data = load_ra_slam_problem(data_path);

    [X, Fval_base_dim, X_is_optimal, Xround, Fval_lifted, manopt_info, Manopt_opts] = ra_slam(exp_data, manopt_opts, do_not_lift);

    res = struct();

    % values from the true final solution
    res.X = X; % the final solution (rounded and refined)
    res.Fval_final = Fval_base_dim; % cost of the final soln
    res.X_is_certified = X_is_optimal;

    % values from rounding the lifted solution
    res.Xround = Xround; % the rounded solution

    % values from the last lifted optimization
    res.Fval_lifted = Fval_lifted; % the final cost of the lifted optimization
    res.time = manopt_info.time; % time elapsed at each iteration

    save_experiment_results(res, manopt_info, data_path);
    plot_solution(X, exp_data)
    animate_lifted_solver_trajectory(data_path);
end

% make a noise when done
load handel
sound(y(1:20000, :),Fs)

