function [least_eigvec, least_eigval] = get_saddle_escape_direction(S, X)
    % using strategy from Rosen, "Accelerating Certifiable Estimation with
    % Preconditioned Eigensolvers"
    assert (size(S,1) == size(S,2), 'S must be square, but is %d x %d', size(S,1), size(S,2));
    n = size(S,1);
    if nargin < 2
        X = rand(n, 5);
    end

    assert (size(X,1) == n, 'X must have same number of rows as S, but is %d x %d', size(X,1), size(X,2));

    x_cols = size(X, 2);
    num_xcols_use = max(2, x_cols-3);
    n_extra_cols = 4;
    block_size = x_cols + n_extra_cols;
    initial_block = zeros(n, block_size);
    initial_block(:, 1:num_xcols_use) = X(:, 1:num_xcols_use);
    initial_block(:, num_xcols_use+1:end) = rand(n, block_size - num_xcols_use);

    global L
    global LT
    global D
    global P
    global PT
    % global p
    % global scale
    [L, D_orig, P] = ldl(S);
    % [L, D_orig, p, scale, ~] = ildl(S);
    LT = L';
    PT = P';
    % Ainv = scale P Linv' Dinv Linv P' scale


    % for each block 'k' of D_orig, let the corresponding block of D be
    % D_k = Q * diag(1/ abs(lambda_1), ..., 1/ abs(lambda_block_size)) * Q'
    % where Q is the matrix of eigenvectors of the block
    % and lambda_i is the i-th eigenvalue of the block.
    % Importantly, D_k is either of block_size 1x1 or 2x2.
    D = zeros(size(D_orig));

    % iterate over the diagonal of D_orig
    cur_idx = 1;
    num_blk_size_one = 0;
    num_blk_size_two = 0;
    while cur_idx <= n
        % if the adjacent off-diagonal entry is zero, then the current
        % diagonal entry is a 1x1 block
        if cur_idx == n || D_orig(cur_idx + 1, cur_idx) == 0
            val = D_orig(cur_idx, cur_idx);
            if val < 0
                abs_val_recip = -1 / val;
            else
                abs_val_recip = 1 / val;
            end
            D(cur_idx, cur_idx) = abs_val_recip;
            cur_idx = cur_idx + 1;
            num_blk_size_one = num_blk_size_one + 1;
        else
            blk = D_orig(cur_idx:cur_idx+1, cur_idx:cur_idx+1);
            [Q, lambda] = eigs(blk);
            D(cur_idx:cur_idx+1, cur_idx:cur_idx+1) = ...
                Q * diag(1 ./ abs(diag(lambda))) * Q';
            cur_idx = cur_idx + 2;
            num_blk_size_two = num_blk_size_two + 1;
        end
    end
    verbosity_level = 0;
    [block_eigvecs, eigvals, failureFlag, lambdaHist, resNormHist] = lobpcg(...
        initial_block, S,...
        [], 'precfun',...
        1e-4, 50,  verbosity_level);
    least_eigvec = block_eigvecs(:, 1);
    least_eigval = eigvals(1);

end