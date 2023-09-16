function V = construct_V_matrix(measurements, var_idx_mapping)
%function V = construct_V_matrix(measurements)
%
% This function constructs and returns the translational data matrix V
% defined in equation (16) of the paper

% Copyright (C) 2016 by David M. Rosen

D = length(measurements.t{1}); % D = dimension of SE(d)
N = var_idx_mapping.pose_var_name_to_cnt.length;  % N = number of nodes in the pose graph
M = size(measurements.edges,1); % M = number of edges in the pose graph

% The number of nonzero elements in V; there are D nonzero elements for
% each translational measurement, plus D*N slots to store the sum of the
% observations emanating from each node

NNZ = D*M + D*N;

rows = zeros(1, NNZ);
cols = zeros(1, NNZ);
vals = zeros(1, NNZ);

out_nodes = cell2mat(values(var_idx_mapping.pose_var_name_to_cnt, measurements.edges(:, 1)));
in_nodes = cell2mat(values(var_idx_mapping.pose_var_name_to_cnt, measurements.edges(:, 2)));

for e = 1:M

    %Extract measurement data
    i = out_nodes(e);  %The node that this edge *leaves*
    j = in_nodes(e);   %The node that this edge *enters*

    tij = measurements.t{e};
    tau_ij = measurements.tau{e};

    %Set V_ji = -tau_ij * t_ij'

    cmin = D*(e-1) + 1;
    cmax = D*(e-1) + D;

    rows(cmin:cmax) = j*ones(1, D);
    cols(cmin:cmax) = D*(i-1) + 1 : D*(i-1) + D;
    vals(cmin:cmax) = -tau_ij * tij';

    %Add this observation to the weighted sum of measurements emanating
    %from node i

    cmin = D*M + D*(i - 1) + 1;
    cmax = D*M + D*(i - 1) + D;
    vals(cmin : cmax) = vals(cmin : cmax) + tau_ij*tij';
end

% Fill in the indices for the running sums

for i = 1:N
    cmin = D*M + D*(i-1) + 1;
    cmax = D*M + D*(i-1) + D;

    rows(cmin:cmax) = i*ones(1,D);
    cols(cmin:cmax) = [D*(i-1) + 1 : D*(i - 1) + D];
end

V = sparse(rows, cols, vals, N, D*N);

