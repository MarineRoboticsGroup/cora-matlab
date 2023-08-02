function [L, LT] = get_ichol_preconditioner(Q)
    ichol_opts.diagcomp = 2e-3;
    L = ichol(Q, ichol_opts);
    LT = L';
end