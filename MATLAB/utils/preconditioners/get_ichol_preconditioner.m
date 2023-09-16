function precon_function = get_ichol_preconditioner(Q)
    ichol_opts.diagcomp = 1e-1;
    L = ichol(Q, ichol_opts);
    LT = L';
    precon_function = @(x) LT \ (L \ x);
end