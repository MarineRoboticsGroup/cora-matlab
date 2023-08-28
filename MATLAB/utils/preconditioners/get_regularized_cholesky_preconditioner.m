function precon_function = get_regularized_cholesky_preconditioner(Q, condition_number_ub)
    [~, neg_greatest_eigval] = get_saddle_escape_direction(-Q);
    if (neg_greatest_eigval > 0)
        error('Failed to compute the greatest eigenvalue of the data matrix Q, should be negative')
    end
    lambda_reg = -neg_greatest_eigval / (condition_number_ub - 1);

    % https://www.mathworks.com/help/matlab/ref/ichol.html#mw_b1e10f55-a0c7-4d70-aedf-32ce0761ac96
    Q_width = size(Q, 1);
    ichol_opts.diagcomp = lambda_reg;
    L = ichol(Q(1:Q_width-1, 1:Q_width-1), ichol_opts);
    LT = L';
    precon_function = @(x) [LT \ (L \ x(1:Q_width-1, :)); zeros(1, size(x, 2))];
end