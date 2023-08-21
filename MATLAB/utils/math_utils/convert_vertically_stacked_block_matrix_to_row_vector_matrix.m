function V = convert_vertically_stacked_block_matrix_to_row_vector_matrix(stacked_block_mat, block_height, blocks_symmetric)
    %{
    Takes a matrix of stacked blocks and converts it such that each column is
    the row-wise vectorized version of the corresponding block

    E.g.,   A = [A1; A2; A3; A4] -> V = [vec(A1), vec(A2), vec(A3), vec(A4)]
            A \in R^{4n x m}
            V \in R^{4n*m x 4}

    If the blocks are symmetric, then we can just take the upper triangular
    part of the matrix and vectorize that
    %}
    if nargin < 3
        blocks_symmetric = false;
    end
    if blocks_symmetric
        % warning("Assuming that the blocks are symmetric");
    end

    stacked_height = size(stacked_block_mat, 1);
    width = size(stacked_block_mat, 2);

    % make sure that the stacked block matrix is a multiple of the block height
    assert(mod(stacked_height, block_height) == 0);

    % get the number of blocks
    num_blocks = stacked_height / block_height;
    if blocks_symmetric
        assert (block_height == width);
        vector_height = (block_height * (block_height+1)) / 2;
    else
        vector_height = block_height * width;
    end

    if issparse(stacked_block_mat)
        % get the indices and values of the matrix
        [row_indices, block_col_indices, values] = find(stacked_block_mat);

        % block row indices = mod(row_indices, block_height)
        block_row_indices = mod(row_indices-1, block_height)+1; % get the row *within* the block

        % if symmetric, we want to ignore any elements where the block_row_index is greater than the block_col_index
        if blocks_symmetric
            % get the indices of the elements that we want to keep
            keep_indices = block_row_indices <= block_col_indices;

            % keep only the elements that we want
            block_row_indices = block_row_indices(keep_indices);
            block_col_indices = block_col_indices(keep_indices);
            values = values(keep_indices);
            row_indices = row_indices(keep_indices);

            % for any values on the block diagonal, we need to divide by 2
            diagonal_indices = block_row_indices == block_col_indices;
            values(diagonal_indices) = values(diagonal_indices) / sqrt(2);
        end

        % calculate the row offset (i.e., how much the vectorized index is
        % offset due to the row stacking)
        if blocks_symmetric
            % if symmetric then the column offset is sum of the lengths of the
            % upper-triangular sections of all previous columns [c_off = c*(c-1)/2]
            row_offset = (block_row_indices-1) .* block_row_indices / 2;
        else
            row_offset = (block_row_indices-1) * width; % get the row offset due to stacking the columns.
        end
        new_row_indices = block_col_indices + row_offset;

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
        error("Dense matrices not supported");
    end


end

