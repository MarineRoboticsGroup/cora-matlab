function precon_function = get_regularized_cholesky_preconditioner(prob, condition_number_ub)
Dnorm = eigs(prob.Q, 1, 'largestreal');
lambda_reg = Dnorm / (condition_number_ub - 1);


if prob.use_marginalized
    idxs = [prob.all_R_idxs, prob.all_d_idxs];
    [R, ~, P] = chol(prob.Q(idxs, idxs) + speye(length(idxs))*lambda_reg, 'lower');
    RT = R';
    PT = P';
    precon_function = @(x) P * (RT \ (R \ (PT * x)));
else
    Q_width = size(prob.Q, 1);
    I_comp = lambda_reg * speye(Q_width-1);
    [R, ~, P] = chol(prob.Q(1:Q_width-1, 1:Q_width-1) + I_comp, 'lower');
    RT = R';
    PT = P';
    
    % Define the preconditioning function
    precon_function = @(x) [P * (RT \ (R \ (PT * x(1:Q_width-1, :))));
        zeros(1, size(x, 2))];
end
end