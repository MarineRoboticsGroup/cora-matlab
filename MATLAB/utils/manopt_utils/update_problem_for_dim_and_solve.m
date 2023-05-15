function [Xlift, Fval, manopt_info, Manopt_opts] = update_problem_for_dim_and_solve(problem, lifted_dim, init_point, Manopt_opts, perturb_lifted_init)
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

    % get the initialization
    function lifted_init = lift_init_point(X, lifted_manifold, add_noise)

        % require that 3 args are given
        if nargin < 3
            error('lift_init_point requires 3 arguments: X, lifted_manifold, add_noise');
        end

        base_dim = lifted_manifold.base_dim;
        lift_dim = lifted_manifold.lifted_dim;

        % make sure base_dim <= size(X, 2) <= lift_dim
        assert(base_dim <= size(X, 2) && size(X, 2) <= lift_dim, ...
            'lift_init_point: X must have size(X, 2) in [base_dim, lifted_dim] but has shape %s', mat2str(size(X)));

        % pad X with zeros to make it the right shape
        lifted_init = lifted_manifold.zeros();
        lifted_init(:, 1:size(X, 2)) = X;

        % randomly apply a rotation to the lifted point to fill the zero entries
        % but keeping the rotation the same (up to gauge symmetry)
        rot = randrot(lift_dim);
        lifted_init = lifted_init * rot;

        % if desired, add some noise to the lifted point. This is useful becausei
        % if rank(lifted_init) < lifted_dim, then it is likely (due to the
        % fundamentals of manifold optimization) that we will not be able to
        % estimate a higher rank point
        if add_noise
            lifted_init = lifted_init + 1e-3 * rand(size(lifted_init));
        end
    end

    % if init_point is [], then get a random point from M
    add_noise = perturb_lifted_init;
    if isempty(init_point)
        X0 = M.rand();
    else
        X0 = lift_init_point(init_point, M, add_noise);
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
