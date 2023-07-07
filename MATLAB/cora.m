function [X, final_soln_optimal, cora_iterates_info, Manopt_opts] = cora(problem, Manopt_opts)

    % assert that 2 args are given
    if nargin < 2
        error('requires 2 arguments: problem and Manopt_opts');
    end


    %%%%%% general problem info
    Q = problem.Q;
    function [f, store] = cost(X, store)
        if ~isfield(store, 'QX')
            store.QX = Q*X;
        end
        QX = store.QX;

        f = 0.5* trace(X'*QX);
    end

    function [g, store] = egrad(X, store)
        if ~isfield(store, 'QX')
            store.QX = Q*X;
        end
        g = store.QX;
    end

    function [h, store] = ehess(~, dX, store)
        h = Q * dX;
    end

    problem.cost = @cost;
    problem.egrad = @egrad;
    problem.ehess = @ehess;

    %%%%%% problem info - unique to lifted dimension

    base_dim = problem.dim;

    % check Manopt_opts.init to see if it is "odom", or "gt"
    if Manopt_opts.init == "odom"
        init_point = problem.X_odom;
    elseif Manopt_opts.init == "gt"
        init_point = problem.X_gt;
    else
        init_point = [];
        warning("No valid initialization given. Using random point as default.")
    end

    % if we know there are no loop closures and we're randomly initializing then
    % let's skip a few dimensions higher (from experience)
    lifted_dim = base_dim; % start out by lifting some dimensions

    cora_iterates_info = [];
    soln_is_optimal = false;
    add_noise_to_lifted_pt = true;
    min_eigvec = [];
    min_eigval = [];
    while ~soln_is_optimal
        % solve the lifted problem and try to certify it
        fprintf("Trying to solve at rank %d\n", lifted_dim);
        [Xlift, Fval_lifted, manopt_info, Manopt_opts] = update_problem_for_dim_and_solve(problem, lifted_dim, init_point, Manopt_opts, add_noise_to_lifted_pt, min_eigvec, min_eigval);
        [soln_is_optimal, ~, min_eigvec, min_eigval] = certify_solution(problem, Xlift, Manopt_opts.verbosity);
        add_noise_to_lifted_pt = false;

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

    fprintf("The staircase algorithm has found an optimal solution with dimension %d.\n", lifted_dim);
    fprintf("Singular values of lifted solution are: %s \n", mat2str(svd(Xlift)));

    % print the rank, singular values, and cost of the solution
    if Manopt_opts.verbosity > 0
        fprintf("Lifted solution has rank %d\n", rank(Xlift));
        fprintf("Cost of lifted solution is %f\n", Fval_lifted);
    end

    Xround = round_solution(Xlift, problem, Manopt_opts.verbosity);

    % refine the rounded solution with one last optimization
    add_noise_to_pt = false;
    [X, Fval_base_dim, soln_manopt_info, ~] = update_problem_for_dim_and_solve(problem, base_dim, Xround', Manopt_opts, add_noise_to_pt, [], []);
    cora_iterates_info = [cora_iterates_info, soln_manopt_info];

    % print the rank, singular values, and cost of the solution
    if Manopt_opts.verbosity > 0
        fprintf("Refined solution has rank %d\n", rank(X));
        fprintf("Cost of refined solution is %f\n", Fval_base_dim);
    end

    gap = Fval_base_dim - Fval_lifted;
    rel_suboptimality = gap / Fval_base_dim;
    if gap < 1e-5
        warning("Final solution is optimal");
        final_soln_optimal = true;
    else
        % print the gap between the final solution and the optimal solution
        warning("Gap between final solution and optimal solution is %f", gap);
        warning("Relative suboptimality is %f%s", rel_suboptimality*100, "%");
        final_soln_optimal = false;
    end

end