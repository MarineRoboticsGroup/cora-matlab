function V = convert_vertically_stacked_block_matrix_to_column_vector_matrix(stacked_block_mat, block_height, blocks_symmetric)
    %{
    Takes a matrix of stacked blocks and converts it such that each column is
    the column-wise vectorized version of the corresponding block

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
    column_vector_shape = [vector_height, num_blocks];

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

        % calculate the column offset (i.e., how much the vectorized index is offset due to the column stacking)
        if blocks_symmetric
            % if symmetric then the column offset is sum of the lengths of the
            % upper-triangular sections of all previous columns [c_off = c*(c-1)/2]
            column_offset = (block_col_indices-1) .* block_col_indices / 2;
        else
            column_offset = (block_col_indices-1) * block_height; % get the row offset due to stacking the columns.
        end
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

        if blocks_symmetric && num_blocks == 1
            % divide diagonal elements by 2
            idxs = 1:block_height;
            diag_idxs = sub2ind(size(V), idxs, idxs);
            V(diag_idxs) = V(diag_idxs) / sqrt(2);

            % get just upper triangular part
            triu_idxs = triu(true(size(V)));
            V = V(triu_idxs);
        else
            error("Dense symmetric matrices not fully supported");
        end
    end


end

