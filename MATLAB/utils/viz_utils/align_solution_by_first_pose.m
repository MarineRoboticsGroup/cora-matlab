function X = align_solution_by_first_pose(X, dim, first_rot_idxs, t_idxs, l_idxs)
    % make sure we have 4 args
    assert (nargin == 5, "align_solution_by_first_pose requires 5 arguments");
    assert (dim == length(first_rot_idxs), "first_rot_idxs should be a vector of indices");

    [nrows, ncols] = size(X);

    % matrix should be short and wide
    assert(nrows < ncols, "There should be more columns than rows, try transposing the matrix");

    % if there is a prior then the first pose is an "origin" frame and can
    % be used to align the solution
    R_origin = X(:, first_rot_idxs);
    t_origin = X(:, t_idxs(1));
    R_inv_origin = R_origin';
    t_inv_origin = - (R_inv_origin * t_origin);
    X = R_inv_origin * X;
    X(:, t_idxs) = X(:, t_idxs) + t_inv_origin;
    X(:, l_idxs) = X(:, l_idxs) + t_inv_origin;

end
