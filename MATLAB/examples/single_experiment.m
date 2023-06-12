clc; clear all; close all;

% the data directory is two levels up from this file, let's get its path
data_dir = fileparts(fileparts(mfilename('fullpath')));

% the experiment we want to run is in <data_dir>/MATLAB/marine/factor_graph.mat
data_path = fullfile(data_dir, 'MATLAB', 'marine', 'factor_graph.mat');
exp_data = load_ra_slam_problem(data_path);

% set the options for the solver
manopt_opts.init = "random";
manopt_opts.verbosity = 2;
manopt_opts.debug = 0;
manopt_opts = get_manopt_opts(manopt_opts);

% run the solver
[X, X_is_optimal, cora_iterates_info, Manopt_opts] = cora(exp_data, manopt_opts);

% values from the true final solution
res = struct();
res.X = X; % the final solution (rounded and refined)
res.X_is_certified = X_is_optimal;

% save the results and write to .tum format for comparison
save_experiment_results(res, cora_iterates_info, data_path);

% visualize the solution
animation_show_gt = true; % show the ground truth trajectory in the animation
animate_lifted_solver_trajectory(data_path, animation_show_gt);
