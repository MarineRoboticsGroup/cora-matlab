function [L, LT] = get_regularized_cholesky_preconditioner(Q, condition_number_ub)
    [~, least_eigval, not_converged] = get_saddle_escape_direction(-Q);
    if not_converged
        error('Failed to compute the smallest eigenvalue of the data matrix Q')
    end
    lambda_reg = -least_eigval / (condition_number_ub - 1);

    % https://www.mathworks.com/help/matlab/ref/ichol.html#mw_b1e10f55-a0c7-4d70-aedf-32ce0761ac96
    ichol_opts.diagcomp = lambda_reg;
    L = ichol(Q, ichol_opts);
    LT = L';
end