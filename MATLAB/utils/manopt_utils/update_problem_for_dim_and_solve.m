function [Xlift, Fval, manopt_info, Manopt_opts] = update_problem_for_dim_and_solve(problem, lifted_dim, init_point, Manopt_opts, add_noise, perturbation, second_order_descent_val)
    % updates the problem to have the given lifted dimension and solves it
    % using manopt
    %
    % Returns:
    %  Xlift: the solution to the problem
    %  Fval: the cost of the solution
    %  manopt_info: the info struct returned by manopt
    %  Manopt_opts: the options struct used by manopt


    M = RaSlamManifoldFactory(problem, lifted_dim);
    problem.M = M;

    if Manopt_opts.verbosity > 0
        fprintf("Manifold typical distance: %f\n", M.typicaldist());
    end

    % use incomplete Cholesky preconditioner
    function [u, store] = preconditioner(X, U, store, ~, mani, L, LT)
        u = LT \ (L \ U);
        u = mani.tangent(X, u);
    end

    LGrho = problem.Q;
    ichol_opts.diagcomp = 1e-3;
    LGrho = LGrho + ichol_opts.diagcomp * speye(size(LGrho));
    L = ichol(LGrho, ichol_opts);
    LT = L';
    problem.precon = @(X, U, store) preconditioner(X, U, store, problem, M, L, LT);

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
