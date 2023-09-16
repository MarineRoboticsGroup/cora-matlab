function X_new = get_lower_rank_solution(X, K)
    % make sure X is tall and skinny
    assert(size(X, 1) > size(X, 2));

    X_new = X;

    sparse_vmat_flag = false;
    SpRight = get_sp_right(K, X_new);
    [SpRightQ, ~, SpRightNullIdxs, ~] = SpRight{:};
    cnt = 0;
    while(~isempty(SpRightNullIdxs))

        % randomly choose idx from SpRightNullIdxs
        idx = randi(length(SpRightNullIdxs));
        % idx = 1;
        v = SpRightQ(:, SpRightNullIdxs(idx));
        vmat = vector_to_symmetric_matrix(v, sparse_vmat_flag);

        % get step size and apply
        max_eigval = eigs(vmat, 1, "largestreal");
        alpha = -1 / max_eigval;
        I = eye(size(vmat));
        delta = I + alpha*vmat;
        [eigvecs_delta, eigvals_delta] = eig(delta);

        % check that all eigvals are non-negative
        eigvals_delta = diag(eigvals_delta);
        assert(all(eigvals_delta >= -1e-6));

        % trim eigvecs_delta and eigvals_delta to ignore zero eigvals
        nonzero_idxs = find(eigvals_delta > 1e-6);
        eigvals_delta = eigvals_delta(nonzero_idxs);
        eigvecs_delta = eigvecs_delta(:, nonzero_idxs);

        % perform square root update of X
        X_new = X_new * eigvecs_delta * diag(sqrt(eigvals_delta));

        % check if constraints are satisfied and cost unchanged
        % new_constraints_val = compute_constraint_vals(X_new, stacked_constraints_vec);
        % new_cost_val = compute_constraint_vals(X_new, Qvec);
        % constraint_diff = new_constraints_val - original_constraint_val;
        % cost_diff = new_cost_val - original_cost_val;
        % constraint_diff_norm = norm(constraint_diff);
        % cost_diff_norm = norm(cost_diff);
        % fprintf("Constraint diff norm: %f, cost diff norm: %f \n", constraint_diff_norm, cost_diff_norm);

        SpRight = get_sp_right(K, X_new);
        [SpRightQ, ~, SpRightNullIdxs, ~] = SpRight{:};

        % % save Xnew for each iteration
        % cnt = cnt + 1;
        % X_new_path = sprintf("X_new_%d.mat", cnt);
        % save(X_new_path, "X_new");

    end
end


function [eigenvector_sensitivities] = compute_eigenvector_sensitivities(K, X)
    Xsp = sparse(X);
    % K = (l+1)*n x n, where l is the number of constraints and n is the number of variables
    % U = n x d
    % KU = [Q*U; A0*U; A1*U; ...; An*U] \in R^{(l+1)*n x d}
    KX = (K * Xsp);
    n = size(X, 1);
    d = size(X, 2);

    % sensitivities for each matrix (Q + constraints) w.r.t. eigenvectors (columns of U)
    sensitivities_height = size(K, 1) / n;
    sensitivities_width = d*(d+1)/2;
    eigenvector_sensitivities = zeros(sensitivities_height, sensitivities_width);

    for block_idx=1:sensitivities_height
        block = (Xsp') * KX((block_idx-1)*n+1:block_idx*n, :);
        block_U = convert_vertically_stacked_block_matrix_to_column_vector_matrix(block, d, true);
        eigenvector_sensitivities(block_idx, :) = block_U(:);
    end
end

function [SpRight] = get_sp_right(K, X)
    eigenvector_sensitivities = compute_eigenvector_sensitivities(K, X);
    % [SpRight, ~] = spspaces(eigenvector_sensitivities', 1);
    [~, SpRight] = spspaces(eigenvector_sensitivities, 2);
end

function [constraint_vals] = compute_constraint_vals(X, stacked_constraints_vec)

    % make sure sizes match
    Xsdp_d = size(X, 1);
    Xsdp_dim = (Xsdp_d+1)*Xsdp_d/2;
    assert(size(stacked_constraints_vec, 1) == Xsdp_dim);

    Xsdp_vec = convert_vertically_stacked_block_matrix_to_column_vector_matrix(sparse(X * X'), size(X, 1), true);
    constraint_vals = 2*(stacked_constraints_vec' * Xsdp_vec);
end