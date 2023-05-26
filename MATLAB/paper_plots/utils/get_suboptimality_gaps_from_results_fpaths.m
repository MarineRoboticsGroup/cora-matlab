function [abs_gaps, rel_gaps] = get_suboptimality_gaps_from_results_fpaths(results_fpaths_pairs)

    [num_pairs, ~] = size(results_fpaths_pairs);
    abs_gaps = zeros(num_pairs, 1);
    rel_gaps = zeros(num_pairs, 1);

    for i = 1:num_pairs
        results_fpath = results_fpaths_pairs{i, 1};
        results = load(results_fpath).results;
        Fval_final = results.Fval_final;
        Fval_lifted = results.Fval_lifted;
        abs_gap = Fval_final - Fval_lifted;
        rel_gap = abs_gap / Fval_lifted;
        abs_gaps(i) = abs_gap;
        rel_gaps(i) = rel_gap;
    end

end