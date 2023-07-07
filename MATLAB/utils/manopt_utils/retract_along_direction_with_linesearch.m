% inspired partially by the SE-sync escape_saddle function
function [new_pt, linesearch_success] = retract_along_direction_with_linesearch(problem, tolgradnorm, pt, search_dir)

    % search_dir and pt must have same number of rows
    [nrows_pt, ncols_pt] = size(pt);
    [nrows_search_dir, ncols_search_dir] = size(search_dir);
    assert(nrows_pt == nrows_search_dir);

    % if ncols_search_dir is one, left pad it with zero columns
    if ncols_search_dir == 1 && ncols_pt > 1
        search_dir = [zeros(nrows_pt, ncols_pt - 1), search_dir];
    end

    % set the minimum and initial stepsize
    typical_dist = problem.M.typicaldist();
    alpha_min = typical_dist * 1e-7;
    alpha = typical_dist * 1e-2;

    % Vectors of trial stepsizes and corresponding function values
    alphas = [];
    fvals = [];
    cost = getCost(problem, pt);

    % Backtracking line search
    while (alpha >= alpha_min)

        % Retract along the given tangent vector using the given stepsize
        ytest = problem.M.retr(pt, search_dir, alpha);

        % Ensure that the trial point Ytest has a lower function value than
        % the current iterate Y, and that the gradient at Ytest is
        % sufficiently large that we will not automatically trigger the
        % gradient tolerance stopping criterion at the next iteration
        fytest = getCost(problem, ytest);
        % check the cost decrease
        rel_cost_decrease = (cost - fytest) / cost;
        if rel_cost_decrease > 1e-4

            % check the gradient norm
            grad_fytest = getGradient(problem, ytest);
            grad_fytest_norm = norm(grad_fytest);
            if grad_fytest_norm > 3*tolgradnorm

                % check the preconditioned gradient norm
                preconditioned_grad_fytest = problem.precon(ytest, grad_fytest, []);
                preconditioned_grad_fytest_norm = norm(preconditioned_grad_fytest);
                if preconditioned_grad_fytest_norm > 3*tolgradnorm
                    % Accept this trial point and return success
                    new_pt = ytest;
                    linesearch_success = true;
                    return;
                end
            end

        end

        % Record trial stepsize and function value
        alphas = [alphas, alpha];
        fvals = [fvals, fytest];
        alpha = alpha / 1.5;

    end

    % If control reaches here, we failed to find a trial point that satisfied
    % *both* the function decrease *and* gradient bounds.  In order to make
    % forward progress, we will fall back to accepting the trial point that
    % simply minimized the objective value, provided that it strictly *decreased*
    % the objective from the current (saddle) point

    % Find minimum function value from among the trial points
    fmin_iter_idx = find(fvals == min(fvals));

    fmin = fvals(fmin_iter_idx);
    amin = alphas(fmin_iter_idx);

    if (fmin < cost)
        % If this trial point strictly decreased the objective value, accept it and
        % return success
        new_pt = problem.M.retr(pt, search_dir, amin);
        linesearch_success = true;
    else
        % NO trial point decreased the objective value: we were unable to escape
        % the saddle point!
        new_pt = pt;
        linesearch_success = false;
    end
end