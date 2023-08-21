%{
    this file exists to let us call CORA from python via the matlab.engine submodule
    eng.cora_python_interface(
        problem_fpath, animation_show_gt, look_for_cached_soln, nargout=0
    )
%}
function cora_python_interface(problem_fpath, show_animation, animation_show_gt, look_for_cached_soln, solve_marginalized_problem, save_iterates_info)

    if ~exist('save_iterates_info', 'var')
        save_iterates_info = true;
    end

    data_path = problem_fpath;
    assert(endsWith(data_path, ".pyfg"), "data_path must be a .pyfg file");

    res_path = strrep(data_path, ".pyfg", "_results.mat");
    cora_iterates_info_path = strrep(data_path, ".pyfg", '_cora_iterates_info.mat');
    data_dir = fileparts(data_path);
    verbose_parsing = false;

    if look_for_cached_soln && exist(res_path, 'file')
        warning("Results already found for experiment %s ... skipping", data_path);
        if show_animation && exist(cora_iterates_info_path, 'file')
            animate_lifted_solver_trajectory(data_path, animation_show_gt);
        end
        return
    end
    prob = parse_pyfg_text(problem_fpath, solve_marginalized_problem, verbose_parsing);

    % set the options for the solver
    manopt_opts.init = "gt";
    manopt_opts.verbosity = 0;
    manopt_opts.debug = 0;
    manopt_opts = get_manopt_opts(manopt_opts);

    % run the solver
    [X, optimality_info, cora_iterates_info, ~] = cora(prob, manopt_opts);

    % values from the true final solution
    res = struct();
    res.X = X; % the final solution (rounded and refined)
    res.certified_lower_bound = optimality_info.certified_lower_bound;
    res.final_soln_cost = optimality_info.final_soln_cost;

    % save the results and write to .tum format for comparison
    if save_iterates_info
        save_experiment_results(data_dir, res, cora_iterates_info);
    else
        save_experiment_results(data_dir, res);
    end
    write_result_to_tum(X, prob, data_dir)

    % visualize the solution
    if show_animation
        animate_lifted_solver_trajectory(data_path, animation_show_gt);
    end
end