function K = get_constraints_as_stacked_block_matrix(problem)
    %{
    Returns a sparse matrix of stacked constraints, where each constraint
    is a (z x z) block matrix. The constraints are stacked vertically.

    K = [K_1; K_2; ...; K_m] where K_i is a (z x z) block matrix.

    We can similarly write K = [Krot; Krange] where Krot is a series of
    rotation constraints and Krange is a series of range constraints. This is
    how we will actually store the constraints in the matrix.

    It's important to note that our indices are ordered as follows:
        - rotations
        - ranges
        - translations
        - landmarks
    %}

    % extract relevant information from problem
    d = problem.dim;
    num_poses = problem.num_poses;
    num_rotation_constraints_per_pose = d * (d + 1) / 2;
    num_rotation_constraints = num_poses * num_rotation_constraints_per_pose;
    num_range_constraints = problem.num_range_measurements;

    % calculate matrix size
    total_num_constraints = num_rotation_constraints + num_range_constraints;

    if problem.use_marginalized
        block_size = size(problem.Qmain, 1); % is a square matrix
    else
        block_size = size(problem.Q, 1); % is a square matrix
    end

    %{
    the rotation constraints will look like the following (for the 2D case):

    K1 = [1, 0;  K2 = [0, 0;    K3 = [0, 0.5;
          0, 0];       0, 1];         0.5, 0];

    %}

    diag_rows = 1:num_poses*d;
    diag_cols = diag_rows;

    % each subsequent row entry should be offset by block_size (b)
    % such that r1 = r1, r2 = r2+b, r3 = r3+2b, etc.
    diag_row_offsets = block_size * (0:num_poses*d-1);
    diag_rows = diag_rows + diag_row_offsets;
    diag_vals = ones(1, length(diag_rows));

    % offdiags are the non-repeating pairs of 1:d (make sure to include
    % both (i, j) and (j, i)
    base_offdiag_rows = [];
    base_offdiag_cols = [];
    for i = 1:d
        for j = i+1:d
            base_offdiag_rows = [base_offdiag_rows, i, j];
            base_offdiag_cols = [base_offdiag_cols, j, i];
        end
    end
    num_offdiag_entries_per_pose = length(base_offdiag_rows);
    num_offdiag_constraints_per_pose = num_offdiag_entries_per_pose / 2;  % divide by 2 since we have (i, j) and (j, i)
    num_offdiag_constraints = num_poses * num_offdiag_constraints_per_pose;

    % for each pose we want to offset the column and row indices of the corresponding
    % entries by d (the dimension of the rotation matrix).
    % This means applying offsets of d * [0, 0, 1, 1, 2, 2, ...] where the number of
    % repeats is the number of offdiagonal entries per pose.
    pose_grouping_indicator = kron(0:num_poses-1, ones(1, num_offdiag_entries_per_pose));
    offdiag_rows = repmat(base_offdiag_rows, 1, num_poses) + d * pose_grouping_indicator;
    offdiag_cols = repmat(base_offdiag_cols, 1, num_poses) + d * pose_grouping_indicator;
    offdiag_vals = 0.5 * ones(1, length(offdiag_rows));

    % now lets construct groupings for each individual constraint. Each constraint has 2 entries
    % so we want to repeat the grouping indicator twice.
    % E.g., [0, 0, 1, 1, ..., m, m] for m offdiagonal constraints.
    constraint_grouping_indicator = kron(0:num_offdiag_constraints-1, ones(1, 2));

    % We only apply this grouping to offset the rows because blocks are stacked vertically.
    % We will simultaneously offset for all of the diagonal blocks above
    num_diagonal_constraints = length(diag_rows);
    offdiag_rows = offdiag_rows + block_size * (constraint_grouping_indicator + num_diagonal_constraints);

    %%%%% form the stacked matrix for rotation constraints %%%%%
    mat_height = total_num_constraints * block_size;
    mat_width = block_size;

    % Now we do the same thing for the range constraints, which are thankfully
    % just 1s on the diagonals. Each entry is a unique constraint and needs
    % to occupy a separate block
    last_rotation_idx = diag_cols(end);
    range_diag_rows = last_rotation_idx + (1:num_range_constraints);
    range_diag_cols = range_diag_rows;

    % we need to offset each row by (number of rotation constraints * block_size)
    range_diag_rows = range_diag_rows + num_rotation_constraints * block_size;

    % now offset each subsequent row entry by block_size
    range_diag_rows = range_diag_rows + block_size * (0:num_range_constraints-1);

    % form the stacked matrix for range constraints
    range_diag_vals = ones(1, length(range_diag_rows));

    % now bring all of the constraints together
    K = sparse([diag_rows, offdiag_rows, range_diag_rows], ...
               [diag_cols, offdiag_cols, range_diag_cols], ...
               [diag_vals, offdiag_vals, range_diag_vals], ...
               mat_height, mat_width);

end
