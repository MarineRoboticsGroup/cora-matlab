function num_poses = get_num_poses_from_fpath_names(results_fpaths_pairs)
    [num_pairs, ~] = size(results_fpaths_pairs);
    num_poses = zeros(num_pairs, 1);
    for i = 1:num_pairs
        [~, name, ~] = fileparts(results_fpaths_pairs{i, 1});

        % example file:
        % factor_graph_4robots_0.5rangeStddev_9500poses_500ranges_0loopClosures_29997seed_results.mat
        % we want to extract the number of poses - 9500 in this case
        leading_str = "Stddev_";
        trailing_str = "poses_";
        start_idx = strfind(name, leading_str) + strlength(leading_str);
        end_idx = strfind(name, trailing_str) - 1;
        char_name = char(name);
        num_poses(i) = str2double(char_name(start_idx:end_idx));
    end
end