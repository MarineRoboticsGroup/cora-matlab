function Xfull = extract_translations_from_marginalized_solution(X, problem)
    % t* = - (X*)' * Qxy * Qyy^{-1};
    assert(problem.Q_is_marginalized, 'Q must be marginalized to extract translations');
    assert(size(X, 1) == size(problem.Qxy, 1), 'X and Qxy must have the same number of rows');

    translations = - (X' * problem.Qxy) / problem.Qyy;
    Xfull = [X; translations'];
end