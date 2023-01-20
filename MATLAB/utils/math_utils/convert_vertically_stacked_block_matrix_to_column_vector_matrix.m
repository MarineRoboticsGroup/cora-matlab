function V = convert_vertically_stacked_block_matrix_to_column_vector_matrix(stacked_block_mat, block_height)
    stacked_height = size(stacked_block_mat, 1);
    width = size(stacked_block_mat, 2);

    % make sure that the stacked block matrix is a multiple of the block height
    assert(mod(stacked_height, block_height) == 0);

    % get the number of blocks
    num_blocks = stacked_height / block_height;
    vector_height = block_height * width;
    column_vector_shape = [vector_height, num_blocks];

    if issparse(stacked_block_mat)
        % get the indices and values of the matrix
        [row_indices, col_indices, values] = find(stacked_block_mat);

        % rows are the row of the block + the offset due to stacking the
        % columns of the block

        % block row indices = mod(row_indices, block_height)
        block_row_indices = mod(row_indices-1, block_height)+1; % get the row *within* the block
        column_offset = (col_indices-1) * block_height; % get the row offset due to stacking the columns
        new_row_indices = block_row_indices + column_offset;

        % new column indices are just the matrix block that the element is in
        new_col_indices = ceil(row_indices / block_height);

        % check that the new row indices are in the correct range
        assert(all(new_row_indices >= 1));
        assert(all(new_row_indices <= vector_height));
        assert(all(new_col_indices >= 1));
        assert(all(new_col_indices <= num_blocks));

        % create the new sparse matrix
        V = sparse(new_row_indices, new_col_indices, values, vector_height, num_blocks);
    else
        % print a warning so we know we're using a dense matrix
        warning("Reshaping a dense matrix to a column vector matrix");

        % take a matrix of stacked blocks and convert it such that each column is
        % the vectorized version of the corresponding block
        int1 = reshape(stacked_block_mat', [width, block_height, num_blocks]);
        int2 = permute(int1, [2, 1, 3]);
        V = reshape(int2, column_vector_shape);
    end


end

