function X = align_solution_by_first_pose(X, problem_data)
    [nrows, ncols] = size(X);
    assert(nrows < ncols, "There should be more columns than rows, try transposing the matrix");

    % get some useful values and run sanity checks
    dim = problem_data.dim;

    % if there is a prior then the first pose is an "origin" frame and can
    % be used to align the solution
    T_origin = X(:, 1:dim+1);
    R_origin = T_origin(1:dim, 1:dim);
    t_origin = T_origin(:, dim+1);
    R_inv_origin = R_origin';
    t_inv_origin = - (R_inv_origin * t_origin);
    X = R_inv_origin * X;
    X(:, problem_data.all_t_idxs) = X(:, problem_data.all_t_idxs) + t_inv_origin;
    X(:, problem_data.all_l_idxs) = X(:, problem_data.all_l_idxs) + t_inv_origin;

end
