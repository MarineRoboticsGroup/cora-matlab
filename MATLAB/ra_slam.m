function [X, final_soln_optimal, cora_iterates_info, Manopt_opts] = ra_slam(problem, Manopt_opts, do_not_lift)

    % assert that 2 args are given
    if nargin < 2
        error('ra_slam requires 2 arguments: problem and Manopt_opts');
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
%         QX = store.QX;
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
    lifted_dim = base_dim + 20; % start out by lifting some dimensions
    if isfield(problem, 'num_loop_closures')
        if problem.num_loop_closures == 0 && isempty(init_point)
                lifted_dim = base_dim + 3;
        end
    end

    if do_not_lift
        warning("do_not_lift == true, just solving problem at base_dim");
        lifted_dim = base_dim;
    end

    cora_iterates_info = [];
    soln_is_optimal = false;
    while ~soln_is_optimal
        % solve the lifted problem and try to certify it
        warning("Trying to solve at rank %d \n", lifted_dim);
%         check_value_is_valid(problem, init_point);
        perturb_lifted_init = true;
        [Xlift, Fval_lifted, manopt_info, Manopt_opts] = update_problem_for_dim_and_solve(problem, lifted_dim, init_point, Manopt_opts,perturb_lifted_init);
        soln_is_optimal = certify_solution(problem, Xlift);

        % add all of the new Xvals from manopt_info to cora_iterates_info but do not
        % add anything else from manopt_info
        cora_iterates_info = [cora_iterates_info, manopt_info];

        if do_not_lift || lifted_dim > 10
            warning("Exiting without finding optimal solution - either do_not_lift is true or lifted_dim > 10");
            soln_is_optimal = 1;
        end

        % if solution is not optimal, increase lifted dimension and try again from previous solution
        if ~soln_is_optimal
            lifted_dim = lifted_dim + 1;
            init_point = Xlift;
        end
    end

    fprintf("The staircase algorithm has found an optimal solution with dimension %d.\n", lifted_dim);

    % print the rank, singular values, and cost of the solution
    fprintf("Lifted solution has rank %d\n", rank(Xlift));
    fprintf("Singular values of lifted solution are: %s \n", mat2str(svd(Xlift)));
    fprintf("Cost of lifted solution is %f\n", Fval_lifted);

    Xround = round_solution(Xlift, problem);

    % refine the rounded solution with one last optimization
    check_value_is_valid(problem, Xround);
    perturb_lifted_init = false;
    [X, Fval_base_dim, soln_manopt_info, ~] = update_problem_for_dim_and_solve(problem, base_dim, Xround', Manopt_opts, perturb_lifted_init);
    cora_iterates_info = [cora_iterates_info, soln_manopt_info];

    % print the rank, singular values, and cost of the solution
    fprintf("Refined solution has rank %d\n", rank(X));
    fprintf("Singular values of refined solution are: %s \n", mat2str(svd(X)));
    fprintf("Cost of refined solution is %f\n", Fval_base_dim);


    % print if the final solution is optimal
    final_soln_optimal = certify_solution(problem, X);
    fprintf("Final solution is optimal: %d\n", final_soln_optimal);

    % print the gap between the final solution and the optimal solution
    fprintf("Gap between final solution and optimal solution is %f\n", Fval_base_dim - Fval_lifted);

end