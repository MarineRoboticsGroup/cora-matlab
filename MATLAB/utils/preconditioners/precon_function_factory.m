function precon_function = precon_function_factory(problem, opts)
    if problem.use_marginalized
        A = problem.Qmain;
    else
        A = problem.Q;
    end

    if strcmp(opts.type, 'block_jacobi')
        precon_function = get_block_jacobi_preconditioner(A, opts.block_size);
    elseif strcmp(opts.type, 'ichol')
        precon_function = get_ichol_preconditioner(A);
    elseif strcmp(opts.type, 'jacobi')
        precon_function = get_jacobi_preconditioner(A);
    elseif strcmp(opts.type, 'regularized_cholesky')
        precon_function = get_regularized_cholesky_preconditioner(A, opts.condition_number_ub);
    elseif strcmp(opts.type, 'none')
        precon_function = @(x) x;
    else
        supported_types = {'block_jacobi', 'ichol', 'jacobi', 'regularized_cholesky', 'none'};
        error('Unknown preconditioner type: %s \nShould be one of: %s',...
            opts.type, strjoin(supported_types, ', '));
    end
end