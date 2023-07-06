function [new_pt] = add_perturbation_with_linesearch(problem, pt, perturbation, df0)
    if nargin < 4
        df0 = getDirectionalDerivative(problem ,pt, perturbation);
    end

    % perturbation and pt must have same number of rows
    [nrows_pt, ncols_pt] = size(pt);
    [nrows_perturbation, ncols_perturbation] = size(perturbation);
    assert(nrows_pt == nrows_perturbation);

    % if ncols_perturbation is one, left pad it with zero columns
    if ncols_perturbation == 1 && ncols_pt > 1
        perturbation = [zeros(nrows_pt, ncols_pt - 1), perturbation];
    end

    cost = getCost(problem, pt);
    [~, new_pt, ~, ~] = linesearch( ...
                             problem, pt, perturbation, cost,...
                             df0);

% Outputs
%
%  stepsize : norm of the vector retracted to reach newx from x.
%  newx : next iterate suggested by the line-search algorithm, such that
%         the retraction at x of the vector alpha*d reaches newx.
%  newkey : key associated to newx in storedb
%  lsstats : statistics about the line-search procedure
%            (stepsize, number of trials etc).
%
end