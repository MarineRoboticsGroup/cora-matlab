function sym_mat = vector_to_symmetric_matrix(vec, sparse_flag)
    if nargin < 2
        sparse_flag = false;
    end

    % This function converts a vector to a symmetric matrix The vector is
    % assumed to be one of the triangular parts of the matrix with the diagonal
    % elements divided by sqrt(2) s.t. the inner product of two vectors is 1/2
    % the inner product of the two matrices

    vec_len = length(vec);

    % vec_len = d*(d+1)/2 = d^2/2 + d/2
    % d^2 + d - 2*vec_len = 0
    d = (-1 + sqrt(1 + 8*vec_len))/2;
    assert (d * (d+1)/2 == vec_len)

    if sparse_flag
        sym_mat = sparse(d,d);
    else
        sym_mat = zeros(d,d);
    end

    idx = 1;
    [row_idxs, col_idxs] = find(triu(true(d)));
    assert (length(row_idxs) == vec_len)

    % scaling = 1 if offdiagonal, sqrt(2) if diagonal
    vec_scaling = (row_idxs == col_idxs) / sqrt(2) + (row_idxs ~= col_idxs);
    sym_mat(sub2ind(size(sym_mat), row_idxs, col_idxs)) = vec .* vec_scaling;

    sym_mat = sym_mat + sym_mat';
end