function [is_opt, min_eigvec, min_eigval] = certify_solution(problem_data, X, verbose)
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
        return
    end


    % if problem_data.use_marginalized
    %     spX = sparse(extract_translations_from_marginalized_solution(X, problem_data));
    % else
    %     spX = sparse(X);
    % end
    spX = sparse(X);
    QX_height = problem_data.block_size;
    vecQX = convert_vertically_stacked_block_matrix_to_column_vector_matrix(Qproduct(spX, problem_data), QX_height);


    % get the constraint gradients
    stacked_constraint_grads = stacked_constraints * spX;
    block_height = size(spX, 1); % this is the width and height of the blocks of stacked_constraints
    assert(block_height == size(stacked_constraints, 2));


    % solve for lambdas as from: vecKX * lambdas = -vecQX
    K = convert_vertically_stacked_block_matrix_to_column_vector_matrix(stacked_constraint_grads, block_height);
    lambdas = K \ -vecQX;

    % vectorize just the constraints and multiply by lambdas to weight the
    % constraints by their Lagrange multipliers
    vecConstraints = convert_vertically_stacked_block_matrix_to_column_vector_matrix(stacked_constraints, block_height);
    weighted_constraints = vecConstraints * lambdas;

    % unstack the columns of the weighted constraints
    weighted_constraints = reshape(weighted_constraints, block_height, []);

    if problem_data.use_marginalized
        S = weighted_constraints + problem_data.Qmain;
    else
        S = weighted_constraints + problem_data.Q;
    end

    % set is_opt to true if we can cholesky factorize S
    beta = 1e-5;
    is_opt = testPSDofMat(S, beta);
    if is_opt
        if verbose
            fprintf("Solution certified with beta = %d \n", beta);
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
        [min_eigvec, min_eigval] = get_saddle_escape_direction(S);
        if min_eigval > 0
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

