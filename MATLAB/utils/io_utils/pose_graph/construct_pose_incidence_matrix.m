function A = construct_pose_incidence_matrix(measurements, var_idx_mapping)
%function A = construct_pose_incidence_matrix(measurements)
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

N = var_idx_mapping.pose_var_name_to_cnt.length;  % N = number of nodes in the pose graph
M = size(measurements.edges, 1);  % Number of edges in the pose graph

out_nodes = cell2mat(values(var_idx_mapping.pose_var_name_to_cnt, measurements.edges(:, 1)));
in_nodes = cell2mat(values(var_idx_mapping.pose_var_name_to_cnt, measurements.edges(:, 2)));

assert(all(out_nodes >= 1) && all(out_nodes <= N));
assert(all(in_nodes >= 1) && all(in_nodes <= N));

node_indices = [out_nodes, in_nodes];
edge_indices = [ [1:M], [1:M] ];
vals = [-ones(1, M), ones(1, M)];

A = sparse(node_indices, edge_indices, vals, N, M)';

end

