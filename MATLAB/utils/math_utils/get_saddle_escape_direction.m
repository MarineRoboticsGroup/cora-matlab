function [least_eigvec, least_eigval, not_converged] = get_saddle_escape_direction(S)
    % using strategy from Rosen, "Accelerating Certifiable Estimation with
    % Preconditioned Eigensolvers"
    n = size(S,1);
    spI = speye(n);

    block_size = 6;
    initial_block = zeros(n, block_size);
    initial_block(:, 1) = 1e-3 * rand(n, 1);
    initial_block(:, 2:end) = rand(n, block_size - 1);

    global L
    global LT
    global D
    global P
    global PT
    nu = 1e-8;
    M = S + (nu * spI);
    [L, D_orig, P] = ldl(M);
    LT = L';
    PT = P';

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
    [block_eigvecs, eigvals, ~, ~, ~] = lobpcg(...
        initial_block, M,...
        [], 'precfun',...
        1e-4, 50,  verbosity_level);

    least_eigvec = block_eigvecs(:, 1);
    least_eigval = eigvals(1);
    not_converged = least_eigval > 0;

end