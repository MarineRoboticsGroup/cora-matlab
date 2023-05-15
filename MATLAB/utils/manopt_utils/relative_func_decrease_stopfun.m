function [stop, reason] = relative_func_decrease_stopfun(manopt_problem, x, info, last, rel_func_decrease_tol)
%function stop = relative_func_decrease_stopfun(manopt_problem, x, info, last, rel_func_decrease_tol)
%
% This function provides an additional stopping criterion for the Manopt
% library: stop whenever the relative decrease in the function value
% between two successive accepted update steps is less than
% rel_func_decrease_tol
% Copyright (C) 2016 by David M. Rosen
% Modified 2023 by Alan Papalia

%%% MANOPT DOCS - https://www.manopt.org/reference/manopt/core/stoppingcriterion.html
%  function [stop, reason] = stoppingcriterion(problem, x, options, info, last)
%
%  Executes standard stopping criterion checks, based on what is defined in
%  the info(last) stats structure and in the options structure.
%
%  The returned number 'stop' is 0 if none of the stopping criteria
%  triggered, and a (strictly) positive integer otherwise. The integer
%  identifies which criterion triggered:
%   0 : Nothing triggered;
%   1 : Cost tolerance reached;
%   2 : Gradient norm tolerance reached;
%   3 : Max time exceeded;
%   4 : Max iteration count reached;
%   6 : User defined stopfun criterion triggered.
%  The output 'reason' is a string describing the triggered event.

this_iterate_accepted = info(last).accepted;

if (~this_iterate_accepted || last == 1)
    stop = 0;
else
    % This iterate was an accepted update step, so get the index of the
    % most recent previously-accepted update

    accepted_array = [info.accepted];
    previous_accepted_iterate_idx = find(accepted_array(1:last-1), 1, 'last');

    if ~isempty(previous_accepted_iterate_idx)
        % Get the function value at the previous accepted iterate

        previous_val = info(previous_accepted_iterate_idx).cost;
        current_val = info(last).cost;

        rel_change_in_val = (previous_val - current_val) / previous_val;
        if rel_change_in_val < rel_func_decrease_tol
            stop = 6;
            % reason = 'relative decrease in function below tolerance';
            % fprintf('Stopping due to relative decrease in function below tolerance: %e < %e\n', rel_change_in_val, rel_func_decrease_tol);
            reason = sprintf('relative decrease in function below tolerance: %e < %e', rel_change_in_val, rel_func_decrease_tol);
        else
            stop = 0;
        end

    else
        % There have been no previously-accepted update steps
        stop = 0;
    end

end
end


