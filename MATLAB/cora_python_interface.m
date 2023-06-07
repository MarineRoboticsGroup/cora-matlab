%{
    this file exists to let us call CORA from python via the matlab.engine submodule
    eng.cora_python_interface(
        problem_fpath, animation_show_gt, look_for_cached_soln, nargout=0
    )
%}
function cora_python_interface(problem_fpath, animation_show_gt, look_for_cached_soln)

    data_path = problem_fpath;
    res_path = strrep(data_path, ".mat", "_results.mat");
    cora_iterates_info_path = strrep(data_path, '.mat', '_cora_iterates_info.mat');
    data_dir = fileparts(data_path);
    exp_data = load_ra_slam_problem(data_path);

    if look_for_cached_soln && exist(res_path, 'file') && exist(cora_iterates_info_path, 'file')
        warning("Results already found for experiment %s ... skipping", data_path);
        res = load(res_path);
        X = res.results.X;
        write_result_to_tum(X, exp_data, data_dir)
        animate_lifted_solver_trajectory(data_path, animation_show_gt);
        return
    end

    exp_data = load_ra_slam_problem(data_path);

    % set the options for the solver
    manopt_opts.init = "random";
    manopt_opts.verbosity = 2;
    manopt_opts.debug = 0;
    manopt_opts = get_manopt_opts(manopt_opts);

    % run the solver
    [X, X_is_optimal, cora_iterates_info, ~] = cora(exp_data, manopt_opts);

    % values from the true final solution
    res = struct();
    res.X = X; % the final solution (rounded and refined)
    res.X_is_certified = X_is_optimal;

    % save the results and write to .tum format for comparison
    save_experiment_results(res, cora_iterates_info, data_path);

    % visualize the solution
    animate_lifted_solver_trajectory(data_path, animation_show_gt);
end