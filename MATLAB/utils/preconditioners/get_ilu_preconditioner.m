function precon_function = get_ilu_preconditioner(A)
    [L, U] = ilu(A);
    precon_function = @(x) L \ (U \ x);
end