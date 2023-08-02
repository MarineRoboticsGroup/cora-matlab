function problem_data = convert_to_marginalized_problem(problem_data)
    % swaps row/col indices of Q s.t. we can easily marginalize out the translations.
    % Right now the ordering is X = [R0, t0, ..., Rn, tn, l0, ..., lk, d0, ..., dl]
    % we want the ordering to be X = [R0, ..., Rn, d0, ..., dl, l0, ..., lk, t0, ..., tn]
    % the indices are all indicated in problem_data in the following fields:
    %   - all_R_idxs: [1 2 3 5 6 7 9 10 11 13 14 15 17 18 19 21 22 23 25 26 27 29 30 31 33 … ]
    %   - all_d_idxs: [7018 7019 7020 7021 7022 7023 7024 7025 7026 7027 7028 7029 7030 … ]
    %   - all_l_idxs: 7017
    %   - all_t_idxs: [4 8 12 16 20 24 28 32 36 40 44 48 52 56 60 64 68 72 76 80 84 88 92 … ]

    unmarginalized_R_idxs = problem_data.all_R_idxs;
    unmarginalized_d_idxs = problem_data.all_d_idxs;
    unmarginalized_l_idxs = problem_data.all_l_idxs;
    unmarginalized_t_idxs = problem_data.all_t_idxs;
    problem_data.original_R_idxs = unmarginalized_R_idxs;
    problem_data.original_l_idxs = unmarginalized_l_idxs;
    problem_data.original_t_idxs = unmarginalized_t_idxs;

    % save these idxs to the problem_data struct
    % problem_data.unmarginalized_R_idxs = unmarginalized_R_idxs;
    % problem_data.unmarginalized_d_idxs = unmarginalized_d_idxs;
    % problem_data.unmarginalized_l_idxs = unmarginalized_l_idxs;
    % problem_data.unmarginalized_t_idxs = unmarginalized_t_idxs;

    full_length = length(unmarginalized_R_idxs) + length(unmarginalized_d_idxs) + length(unmarginalized_l_idxs) + length(unmarginalized_t_idxs);
    assert (unmarginalized_d_idxs(end) == full_length);

    % overwrite the Q matrix with the new ordering
    reordered_idxs = [unmarginalized_R_idxs, unmarginalized_d_idxs, unmarginalized_l_idxs, unmarginalized_t_idxs];
    problem_data.Q = problem_data.Q(reordered_idxs, reordered_idxs);

    % overwrite the variable indices with the new ordering
    problem_data.all_R_idxs = 1:length(unmarginalized_R_idxs);
    problem_data.all_d_idxs = problem_data.all_R_idxs(end) + (1:length(unmarginalized_d_idxs));
    problem_data.all_l_idxs = problem_data.all_d_idxs(end) + (1:length(unmarginalized_l_idxs));
    problem_data.all_t_idxs = problem_data.all_l_idxs(end) + (1:length(unmarginalized_t_idxs));

    % marginalize out the translations (and landmarks) and drop near-zero vals (truncate)
    [Qxx, Qxy, Qyy] = marginalize_translations(problem_data);

    problem_data.Qxx = Qxx;
    problem_data.Qxy = Qxy;
    problem_data.Qyy = Qyy;

    % indicate that Q is now marginalized
    problem_data.Q_is_marginalized = true;

    % Qmarg = Qxx - Qxy*(Qyy\Qxy');
    % add a function to access Qmarg without explicitly computing it
    problem_data.Qmarg = @(x) (problem_data.Qxx * x) - problem_data.Qxy*(problem_data.Qyy\(problem_data.Qxy' * x));

    % remove {all_l_idxs, all_t_idxs} from problem_data
    problem_data = rmfield(problem_data, 'all_l_idxs');
    problem_data = rmfield(problem_data, 'all_t_idxs');

    % move the fields X_odom and X_gt to X_..._original
    problem_data.X_odom_original = problem_data.X_odom;
    problem_data.X_gt_original = problem_data.X_gt;

    % overwrite X_odom and X_gt with the new ordering (excluding the translations)
    problem_data.X_odom = problem_data.X_odom([unmarginalized_R_idxs, unmarginalized_d_idxs], :);
    problem_data.X_gt = problem_data.X_gt([unmarginalized_R_idxs, unmarginalized_d_idxs], :);

    % move the field stacked_constraints to stacked_constraints_original
    % problem_data.stacked_constraints_original = problem_data.stacked_constraints;

    % stacked_constraints is a little tricky, because it is a block-column matrix (M = [A; B; C; D])
    % and we want to get just the R_idxs and d_idxs for each block.
    % So we need to do some indexing magic to get the right rows and columns.
    original_block_size = size(problem_data.stacked_constraints, 2); % square matrix
    num_blocks = size(problem_data.stacked_constraints, 1) / original_block_size;
    block_row_idxs = reordered_idxs;
    new_block_size= numel(block_row_idxs);
    new_stack_row_idxs = repmat(int64(0:num_blocks-1), new_block_size, 1) * original_block_size;
    new_stack_row_idxs = new_stack_row_idxs + block_row_idxs';
    new_stack_row_idxs = new_stack_row_idxs(:);

    % make sure new_stack_row_idxs has only unique values
    assert (length(new_stack_row_idxs) == length(unique(new_stack_row_idxs)));

    new_stack_col_idxs = block_row_idxs;
    problem_data.stacked_constraints = problem_data.stacked_constraints(new_stack_row_idxs, new_stack_col_idxs);

end

function [Qxx, Qxy, Qyy] = marginalize_translations(problem_data)
    % [ Q_XX   Q_XY ]
    % [ Q_YX   Q_YY ]
    % S = Q_XX - Q_XY * Q_YY^(-1) * Q_YX

    first_trans_idx = problem_data.all_d_idxs(end)+1;

    % get all of the submatrix of Q above the first translation
    Qxy = problem_data.Q(1:first_trans_idx-1, first_trans_idx:end);
    Qxx = problem_data.Q(1:first_trans_idx-1, 1:first_trans_idx-1);
    Qyy = problem_data.Q(first_trans_idx:end, first_trans_idx:end);

    % Qmarg = Qxx - Qxy*(Qyy\Qxy');
    % assert (size(Qmarg, 1) == size(Qmarg, 2));
    % assert (size(Qmarg, 1) == problem_data.all_d_idxs(end));

end