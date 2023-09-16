function Xfull = extract_translations_from_marginalized_solution(X, problem)
    % t* = - (X*)' * Qxy * Qyy^{-1};
    assert(problem.use_marginalized, 'Q must be marginalized to extract translations');

    translations = - (X' * problem.LeftOperator) / problem.Ltrans;
    Xfull = [X; translations'];
end