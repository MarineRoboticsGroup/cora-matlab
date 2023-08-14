function [A, W, D] = construct_range_graph_matrices(measurements, var_idx_mapping)
    % function A = construct_range_incidence_matrix(measurements)
    %
    % This function computes and returns the oriented incidence matrix of the
    % underlying directed graph of measurements:
    %
    % A_{ie} =  -1, if edge e *leaves* node i,
    %           +1, if edge e *enters* node i,
    %           0, otherwise.
    %
    % (see eq. (7) in the paper).

    % Copyright (C) 2016 by David M. Rosen

    M = size(measurements.edges, 1);  % Number of edges in the pose graph
    num_pose_vars = var_idx_mapping.pose_var_name_to_cnt.length;
    num_landmark_vars = var_idx_mapping.landmark_name_to_cnt.length;
    N = num_landmark_vars + num_pose_vars;

    % join the two mappings (poses and landmarks) into one
    % we want the to order the translations to be (pose translations, landmark translations)
    % so we will offset the landmark translations by the number of pose translations

    all_keys = [var_idx_mapping.pose_var_name_to_cnt.keys, var_idx_mapping.landmark_name_to_cnt.keys];
    offset_landmark_vals = cellfun(@(x) x + num_pose_vars, var_idx_mapping.landmark_name_to_cnt.values, 'UniformOutput', false);
    all_values = [var_idx_mapping.pose_var_name_to_cnt.values, offset_landmark_vals];
    joined_mapping = containers.Map(all_keys, all_values);

    out_nodes = cell2mat(values(joined_mapping, measurements.edges(:, 1)));
    in_nodes = cell2mat(values(joined_mapping, measurements.edges(:, 2)));

    % make sure all out_nodes and in_nodes are in the range [1, N]
    assert(all(out_nodes >= 1) && all(out_nodes <= N));
    assert(all(in_nodes >= 1) && all(in_nodes <= N));

    node_indices = [out_nodes, in_nodes];
    edge_indices = [ [1:M], [1:M] ];
    vals = [-ones(1, M), ones(1, M)];

    % the incidence matrix is sparse, so we use the sparse constructor
    A = sparse(node_indices, edge_indices, vals, N, M)';

    % form the sparse diagonal matrices D and W
    D = spdiags(cell2mat(measurements.range)', 0, M, M);
    W = spdiags(cell2mat(measurements.precision)', 0, M, M);

end
