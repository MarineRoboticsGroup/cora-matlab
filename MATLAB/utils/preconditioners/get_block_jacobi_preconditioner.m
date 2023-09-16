function precon_function = get_block_jacobi_preconditioner(A, block_size)

    % Get the size of the matrix
    n = size(A, 1);

    % have a remainder if the block size does not divide the matrix size
    block_size_remainder = mod(n, block_size);

    % Get the number of blocks
    num_default_blocks = (n - block_size_remainder) / block_size;
    block_jacobi = sparse(n, n);

    % for each block, compute the inverse of the block
    for i = 1:num_default_blocks
        % Get the indices of the block
        block_indices = (i-1)*block_size+1:i*block_size;

        % Get the block
        block = A(block_indices, block_indices);

        % Compute the inverse of the block
        block_jacobi(block_indices, block_indices) = inv(block);
    end

    % if there is a remainder, compute the inverse of the remainder block
    if block_size_remainder > 0
        % Get the indices of the block
        block_indices = num_default_blocks*block_size+1:block_size_remainder;

        % Get the block
        block = A(block_indices, block_indices);

        % Compute the inverse of the block
        block_jacobi(block_indices, block_indices) = inv(block);
    end

    % Return the preconditioner function
    precon_function = @(x) block_jacobi * x;

end
