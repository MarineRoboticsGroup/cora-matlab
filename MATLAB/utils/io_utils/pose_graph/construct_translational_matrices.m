function [T, Omega] = construct_translational_matrices(measurements, var_idx_mapping)
%function [T, omega] = construct_translational_matrices(measurements)
%
% This function computes and returns the matrix T containing the raw
% translational measurements (see eq. (24) in the paper) and the diagonal
% matrix Omega containing the translational measurement precisions on its
% main diagonal (see eq. (23) in the paper).

% Copyright (C) 2016 by David M. Rosen

D = length(measurements.t{1}); % D = dimension of SE(d)
N = var_idx_mapping.pose_var_name_to_cnt.length;  % N = number of nodes in the pose graph
M = size(measurements.edges,1); % M = number of edges in the pose graph

%Allocate storate for sparse matrix
rows = zeros(1, D*M);
cols = zeros(1, D*M);
vals = zeros(1, D*M);

omega = zeros(M, 1);
out_nodes = cell2mat(values(var_idx_mapping.pose_var_name_to_cnt, measurements.edges(:, 1)));
assert(all(out_nodes >= 1) && all(out_nodes <= N));

%Iterate over the measurements in the pose graph
for e = 1:M

    %EXTRACT MEASUREMENT DATA
    k = out_nodes(e);  %The index of the node at the tail of this edge

    tij = measurements.t{e};  %The translation corresponding to this observation

    omega(e) = measurements.tau{e};  %The precision for this translational observation


    %PROCESS MEASUREMENT DATA

    rows(D*(e-1) + 1 : D*(e-1) + D) = e*ones(1,D);
    cols(D*(e-1) + 1 : D*(e-1) + D) = [D*(k-1) + 1 : D*(k-1) + D];
    vals(D*(e-1) + 1 : D*(e-1) + D) = -tij';
end

T = sparse(rows, cols, vals, M, D*N);
Omega = spdiags(omega, 0, M, M);
end

