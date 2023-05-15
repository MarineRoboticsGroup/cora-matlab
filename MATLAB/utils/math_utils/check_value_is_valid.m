function is_valid = check_value_is_valid(problem, X)
    % innocent until proven guilty
    is_valid = true;
    
    % if X is tall and skinny, transpose it
    if size(X, 1) > size(X, 2)
        X = X';
    end


    % for each pose, check if the rotation is valid (i.e., R.T * R = I and
    % det(R) = 1)
    all_rot_idxs = problem.all_R_idxs;
    dim = problem.dim;
    num_poses = problem.num_poses;
    for i = 1:num_poses
        rot_i_idx_start = dim * (i-1) + 1;
        rot_i_idx_end = (dim * i);
        rot_start_idx = all_rot_idxs(rot_i_idx_start);
        rot_end_idx = all_rot_idxs(rot_i_idx_end);
        R = X(:, rot_start_idx:rot_end_idx);
        assert (size(R, 1) == dim);
        assert (size(R, 2) == dim);

        det_R = det(R);
        if (abs(det_R - 1) > 1e-6)
            is_valid = false;
            fprintf('rot %d: det(R) is not 1, det(R) = %f\n',i, det_R);
            return;
        end

        R_T = R';
        R_T_R = R_T * R;
        if (norm(R_T_R - eye(dim)) > 1e-6)
            is_valid = false;
            fprintf('rot %d: R.T * R is not I, residual = %f\n',i, norm(R_T_R - eye(dim)));
            disp(R_T_R);
            return;
        end

    end

    % for each distance variable make sure it is a unit vector
    all_dist_idxs = problem.all_d_idxs;
    for i = 1:length(all_dist_idxs)
        dist_i_idx = all_dist_idxs(i);
        dist_i = X(:, dist_i_idx);
        if (norm(dist_i) > 1 + 1e-6)
            is_valid = false;
            fprintf('dist %d: norm is not 1, norm = %f',i, norm(dist_i));
            return;
        end
    end
end