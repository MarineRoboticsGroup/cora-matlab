function Manopt_opts = get_manopt_opts(Manopt_opts)

    if nargin == 0
        disp('Using default settings for Manopt:');
        Manopt_opts = struct;  % Create empty structure
    else
        disp('Manopt settings:');
    end


    if isfield(Manopt_opts, 'verbosity')
        fprintf(' Solver verbosity: %g\n', Manopt_opts.verbosity);
    else
        Manopt_opts.verbosity = 1;
        fprintf(' Setting solver verbosity: %g [default]\n', Manopt_opts.verbosity);
    end

    if isfield(Manopt_opts, 'init')
        fprintf(' Initialization strategy: %s\n', Manopt_opts.init);
    else
        Manopt_opts.init = "odom";
        fprintf(' Setting initialization strategy: %s [default]\n', Manopt_opts.init);
    end

    if isfield(Manopt_opts, 'tolgradnorm')
        fprintf(' Stopping tolerance for norm of Riemannian gradient: %g\n', Manopt_opts.tolgradnorm);
    else
        Manopt_opts.tolgradnorm = 1e-4;
        fprintf(' Setting stopping tolerance for norm of Riemannian gradient to: %g [default]\n', Manopt_opts.tolgradnorm);
    end

    if isfield(Manopt_opts, 'rel_func_tol')
        fprintf(' Stopping tolerance for relative function decrease: %g\n', Manopt_opts.rel_func_tol);
    else
        Manopt_opts.rel_func_tol = 1e-8;
        warning("VERY LOW RELATIVE FUNCTION TOLERANCE")
        fprintf(' Setting stopping tolerance for relative function decrease to: %g [default]\n', Manopt_opts.rel_func_tol);
    end

    if isfield(Manopt_opts, 'maxinner')
        fprintf(' Maximum number of Hessian-vector products to evaluate in each truncated Newton iteration: %d\n', Manopt_opts.maxinner);
    else
        Manopt_opts.maxinner = 1000;
        fprintf(' Setting maximum number of Hessian-vector products to evaluate in each truncated Newton iteration to: %d [default]\n', Manopt_opts.maxinner);
    end

    if isfield(Manopt_opts, 'miniter')
        fprintf(' Minimum number of trust-region iterations: %d\n', Manopt_opts.miniter);
    else
        Manopt_opts.miniter = 1;
        fprintf(' Setting minimum number of trust-region iterations to: %d [default]\n', Manopt_opts.miniter);
    end

    if isfield(Manopt_opts, 'maxiter')
        fprintf(' Maximum number of trust-region iterations: %d\n', Manopt_opts.maxiter);
    else
        Manopt_opts.maxiter = 500;
        fprintf(' Setting maximum number of trust-region iterations to: %d [default]\n', Manopt_opts.maxiter);
    end

    if isfield(Manopt_opts, 'maxtime')
        fprintf(' Maximum permissible elapsed computation time [sec]: %g\n', Manopt_opts.maxtime);
    end

    % Check if a solver was explicitly supplied
    if(~isfield(Manopt_opts, 'solver'))
        % Use the trust-region solver by default
        Manopt_opts.solver = @trustregions;
    end
    solver_name = func2str(Manopt_opts.solver);
    if (~strcmp(solver_name, 'trustregions') && ~strcmp(solver_name, 'conjugategradient') && ~strcmp(solver_name, 'steepestdescent'))
        error(sprintf('Unrecognized Manopt solver: %s', solver_name));
    end
    fprintf('\nSolving Riemannian optimization problems using Manopt''s "%s" solver\n\n', solver_name);


    % Set additional stopping criterion for Manopt: stop if the relative
    % decrease in function value between successive iterates drops below the
    % threshold specified in SE_Sync_opts.relative_func_decrease_tol
    if(strcmp(solver_name, 'trustregions'))
        Manopt_opts.stopfun = @(manopt_problem, x, info, last) relative_func_decrease_stopfun(manopt_problem, x, info, last, Manopt_opts.rel_func_tol);
    end

    % Log the sequence of iterates visited by the Riemannian Staircase
    Manopt_opts.statsfun = @log_iterates;

end