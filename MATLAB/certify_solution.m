function [is_opt, cert_time, min_eigvec, min_eigval] = certify_solution(problem_data, X, verbose)
    if verbose
        tic
    end

    % set up the problem
    stacked_constraints = problem_data.stacked_constraints;
    min_eigvec = [];
    min_eigval = [];
    % if stacked constraints is empty, then we have no constraints and the
    % solution is optimal (we just solved a convex QP)
    if isempty(stacked_constraints)
        is_opt = true;
        cert_time = toc;
        return
    end

    Q = problem_data.Q;
    spX = sparse(X);

    % get the constraint gradients
    stacked_constraint_grads = stacked_constraints * spX;
    block_height = size(X, 1); % this is the width and height of the blocks of stacked_constraints
    assert(block_height == size(stacked_constraints, 2));
    K = convert_vertically_stacked_block_matrix_to_column_vector_matrix(stacked_constraint_grads, block_height);
    QX = Q * spX;
    QX_height = size(QX, 1);
    vecQX = convert_vertically_stacked_block_matrix_to_column_vector_matrix(QX, QX_height);

    % solve for lambdas as from: vecKX * lambdas = -vecQX
    lambdas = K \ -vecQX;

    % vectorize just the constraints and multiply by lambdas to weight the
    % constraints by their Lagrange multipliers
    vecConstraints = convert_vertically_stacked_block_matrix_to_column_vector_matrix(stacked_constraints, block_height);
    weighted_constraints = vecConstraints * lambdas;

    % unstack the columns of the weighted constraints
    weighted_constraints = reshape(weighted_constraints, block_height, []);

    S = weighted_constraints + Q;

    % set is_opt to true if we can cholesky factorize S
    % iterate from 1e-10, 1e-9, ...
    for i = -5:-5
        beta = 10^(i);
        is_opt = testPSDofMat(S, beta);
        if is_opt
            if verbose
                fprintf("Solution certified with beta = %d \n", beta);
            end
            break
        end
    end

    % time for certification
    if verbose
        cert_time = toc;
        fprintf("Certification took %f seconds \n", cert_time);
    end

    if ~is_opt
        if verbose
            fprintf("Not certified after beta of %d \n", beta);
            tic;
        end

        % get minimum eigenvector of S, which is real and symmetric (Hermitian)
        [min_eigvec, min_eigval, not_converged] = get_saddle_escape_direction(S);
        if not_converged
            warning("Saddle escape search did not converge");
            min_eigvec = [];
            min_eigval = [];
        end
        if verbose
            saddle_time = toc;
            fprintf("Saddle escape search took %f seconds \n", saddle_time);
        end
    end

end

function isPSD = testPSDofMat(mat, reg_term)
    try
        chol(mat + reg_term * speye(size(mat)));
        isPSD = true;
    catch
        isPSD = false;
    end

end

