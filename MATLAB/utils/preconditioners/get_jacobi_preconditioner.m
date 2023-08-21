function precon_function = get_jacobi_preconditioner(A)
    diagonal = diag(A);

    % Create a diagonal preconditioning matrix D with diagonal elements as reciprocals
    D = spdiags(1 ./ diagonal, 0, size(A, 1), size(A, 2));

    % Create a function handle to the preconditioner
    precon_function = @(x) D * x;
end