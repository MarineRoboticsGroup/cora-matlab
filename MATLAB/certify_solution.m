function [is_opt, min_eigvec, min_eigval] = certify_solution(problem_data, X, verbose, certEpsilon)
    if nargin < 3
        verbose = false;
    end
    if nargin < 4
        certEpsilon = 1e-5;
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
        S = problem_data.Q;
        S(1:block_height, 1:block_height) = S(1:block_height, 1:block_height) + weighted_constraints;
    else
        S = weighted_constraints + problem_data.Q ;
    end
    S = S + (certEpsilon * speye(size(S)));

    % set is_opt to true if we can cholesky factorize S
    % is_opt = testPSDofMat(S);
    % warning("Disabled saddle escape search");
    % min_eigvec = [];
    % min_eigval = [];
    % is_opt = testPSDofMat(S);

    if ~problem_data.use_marginalized
        [min_eigvec, min_eigval] = get_saddle_escape_direction(S, X);
    else
        [min_eigvec, min_eigval] = get_saddle_escape_direction(S);
    end
    is_opt = min_eigval >= 0;
    is_opt = testPSDofMat(S);
    if is_opt
        fprintf("Solution certified with certEpsilon = %d \n", certEpsilon);
    else
        if verbose
            fprintf("Not certified after certEpsilon of %d \n", certEpsilon);
        end

        % get minimum eigenvector of S, which is real and symmetric (Hermitian)
        % [min_eigvec, min_eigval] = get_saddle_escape_direction(S);
        if min_eigval > 0
            warning("Saddle escape search did not converge");
            min_eigvec = [];
            min_eigval = [];
        elseif problem_data.use_marginalized
            min_eigvec = min_eigvec(1:problem_data.block_size);
        end
    end
end

function isPSD = testPSDofMat(mat)
    try
        chol(mat);
        isPSD = true;
    catch
        isPSD = false;
    end
end

