function [X, optimality_info, cora_iterates_info, Manopt_opts] = cora(problem, Manopt_opts, preconditioner_type)

    % assert that 2 args are given
    if nargin < 2
        error('requires 2 arguments: problem and Manopt_opts');
    end

    % if preconditioner_type is not given, default to "jacobi"
    if nargin < 3
        preconditioner_type = "block_cholesky";
    end



    %%%%%% general problem info
    function [f, store] = cost(X, store)
        if ~isfield(store, 'QX')
            store.QX = Qproduct(X, problem);
        end
        f = 0.5* trace(X'*store.QX);
    end

    function [g, store] = egrad(X, store)
        if ~isfield(store, 'QX')
            store.QX = Qproduct(X, problem);
        end
        g = store.QX;
    end

    function [h, store] = ehess(~, dX, store)
        h = Qproduct(dX, problem);
    end

    problem.cost = @cost;
    problem.egrad = @egrad;
    problem.ehess = @ehess;

    fprintf("Setting up preconditioner...\n")
    precon_opts = struct();
    precon_opts.type = preconditioner_type;
    precon_opts.condition_number_ub = 1e4;
    precon_opts.block_size = problem.dim+2;
    problem.precon_function = precon_function_factory(problem, precon_opts);

    % check Manopt_opts.init to see if it is "odom", or "gt"
    if Manopt_opts.init == "odom"
        init_point = problem.X_odom;
    elseif Manopt_opts.init == "gt"
        init_point = problem.X_gt;
    else
        init_point = [];
        warning("No valid initialization given. Using random point as default.")
    end

    if problem.use_marginalized && ~isempty(init_point)
        marginalized_dim = problem.dim * problem.num_poses + problem.num_range_measurements;
        init_point = init_point(1: marginalized_dim, :);
    end

    % we can start out by lifting some dimensions if we want
    lifted_dim = problem.dim;

    cora_iterates_info = [];
    soln_is_optimal = false;
    add_noise_to_lifted_pt = true;
    min_eigvec = [];
    min_eigval = [];
    while ~soln_is_optimal
        % solve the lifted problem and try to certify it
        fprintf("Trying to solve at rank %d\n", lifted_dim);
        [Xlift, Fval_lifted, manopt_info, Manopt_opts] = update_problem_for_dim_and_solve(problem, lifted_dim, init_point, Manopt_opts, add_noise_to_lifted_pt, min_eigvec, min_eigval);
        certEpsilon = Fval_lifted * 1e-6;
        [soln_is_optimal, min_eigvec, min_eigval] = certify_solution(problem, Xlift, Manopt_opts.verbosity, certEpsilon);

        % Xlift should be tall and skinny
        assert(size(Xlift, 1) > size(Xlift, 2));

        % if eigpairs are empty then add noise to lifted point
        if isempty(min_eigvec) && isempty(min_eigval)
            add_noise_to_lifted_pt = true;
        else
            add_noise_to_lifted_pt = false;
        end

        % add all of the new Xvals from manopt_info to cora_iterates_info but do not
        % add anything else from manopt_info
        cora_iterates_info = [cora_iterates_info, manopt_info];

        max_lifted_dim = 20;
        if lifted_dim > max_lifted_dim
            warning("Exiting without finding optimal solution - lifted_dim > max_lifted_dim");
            soln_is_optimal = 1;
        end

        % if solution is not optimal, increase lifted dimension and try again from previous solution
        if ~soln_is_optimal
            lifted_dim = lifted_dim + 1;
            init_point = Xlift;
        end
    end

    fprintf("The staircase algorithm has found an optimal solution with dimension %d and cost %f.\n", lifted_dim, Fval_lifted);
    fprintf("Singular values of lifted solution are: %s \n", mat2str(svd(Xlift)));

    % print the rank, singular values, and cost of the solution
    certified_lower_bound = Fval_lifted;

    if size(Xlift, 2) == problem.dim
        X = Xlift;
        final_soln_cost = certified_lower_bound;
        fprintf("Problem was solved at the base dimension, no rounding and refinement needed.\n")
    else
        [X, Fval_base_dim, soln_manopt_info] = round_lifted_solution_and_refine(Xlift, problem, Manopt_opts);
        fprintf("Refined cost is %f\n", Fval_base_dim);
        cora_iterates_info = [cora_iterates_info, soln_manopt_info];
        final_soln_cost = Fval_base_dim;
        gap = final_soln_cost - certified_lower_bound;
        rel_suboptimality = gap / certified_lower_bound;
        if abs(rel_suboptimality) < 0.01
            warning("Final solution is optimal");
        elseif gap < 0
            warning("Suboptimality gap is negative - something is wrong. Gap is %f", gap);
        else
            warning("Gap between final solution and optimal solution is %f", gap);
            warning("Relative suboptimality is %f%s", rel_suboptimality*100, "%");
        end

    end

    optimality_info.certified_lower_bound = certified_lower_bound;
    optimality_info.final_soln_cost = final_soln_cost;
end

function [X, Fval_base_dim, soln_manopt_info] = round_lifted_solution_and_refine(Xlift, problem, Manopt_opts)
    Xround = round_solution(Xlift, problem, Manopt_opts.verbosity);
    assert(check_value_is_valid(problem, Xround))

    % refine the rounded solution with one last optimization
    add_noise_to_pt = false;
    [X, Fval_base_dim, soln_manopt_info, ~] = update_problem_for_dim_and_solve(problem, problem.dim, Xround', Manopt_opts, add_noise_to_pt, [], []);
    assert(check_value_is_valid(problem, X))

    % print the rank, singular values, and cost of the solution
    if Manopt_opts.verbosity > 0
        fprintf("Refined solution has rank %d\n", rank(X));
        fprintf("Cost of refined solution is %f\n", Fval_base_dim);
    end

end