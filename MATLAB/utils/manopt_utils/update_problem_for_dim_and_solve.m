function [Xlift, Fval, manopt_info, Manopt_opts] = update_problem_for_dim_and_solve(problem, lifted_dim, init_point, Manopt_opts, add_noise, perturbation, second_order_descent_val)
    % updates the problem to have the given lifted dimension and solves it
    % using manopt
    %
    % Returns:
    %  Xlift: the solution to the problem
    %  Fval: the cost of the solution
    %  manopt_info: the info struct returned by manopt
    %  Manopt_opts: the options struct used by manopt


    if problem.use_marginalized
        M = MarginalizedRaSlamManifoldFactory(problem, lifted_dim);

    else
        M = RaSlamManifoldFactory(problem, lifted_dim);
    end
    problem.M = M;

    % use incomplete Cholesky preconditioner
    function [u, store] = preconditioner(X, U, store, ~, mani, L, LT)
        u = LT \ (L \ U);
        u = mani.tangent(X, u);
    end

    problem.precon = @(X, U, store) preconditioner(X, U, store, problem, M, problem.L, problem.LT);

    % if init_point is [], then get a random point from M
    if isempty(init_point)
        X0 = M.rand();
    else
        X0 = lift_init_point(problem, init_point, M, add_noise, perturbation, second_order_descent_val, Manopt_opts.tolgradnorm);
    end

    % if no Delta_bar is given, then base it on the typical distance of the manifold
    if ~isfield(Manopt_opts, 'Delta_bar')
        Manopt_opts.Delta_bar = M.typicaldist() * 200;
    end

    % print info on problem.M
    if Manopt_opts.verbosity > 0
        fprintf("Refining solution over manifold %s\n", problem.M.name());
    end
    [Xlift, Fval, manopt_info, Manopt_opts] = manoptsolve(problem, X0, Manopt_opts);
end
