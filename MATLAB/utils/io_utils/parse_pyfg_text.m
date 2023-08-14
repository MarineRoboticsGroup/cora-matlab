function problem = parse_pyfg_text(pyfg_fpath, use_marginalized)
    %{
    Parses a pyfg text file into a problem struct.

    Expected example lines:
        VERTEX_SE2 <timestamp> <VertexName> <x> <y> <theta>
        VERTEX_SE2 3152.0 A0 -34.2086489999201 45.300763999112 1.1205036535897932

        VERTEX_XY <timestamp> <VertexName> <x> <y>
        VERTEX_XY L3 1.70946300006472 -5.81220300029963

        EDGE_SE2 <timestamp> <VertexName1> <VertexName2> <x> <y> <theta> <information>
        EDGE_SE2 3511.29851794243 A3587 A3588 0.332397242830102 0.0 -0.0319438761202058 0.010000000000000002 0.0 -0.0 0.010000000000000002 -0.0 0.0001

        EDGE_RANGE <timestamp> <VertexName1> <VertexName2> <range> <stddev>
        EDGE_RANGE 3542.90110708028 A3904 L3 26.69216622156613 0.2655983804002221
    %}
    % make sure the file exists
    assert(exist(pyfg_fpath, 'file') == 2, 'File does not exist: %s', pyfg_fpath);
    assert(nargin == 2, 'Must provide use_marginalized flag as second argument.')

    problem.use_marginalized = use_marginalized;

    if use_marginalized
        fprintf('Using marginalized product.\n');
    else
        fprintf('Not using marginalized product.\n');
    end

    %%%%% Parse the file %%%%%
    tic;
    % [pose_measurements, range_measurements, var_idx_mapping, dim, true_rots, true_pose_translations, true_landmarks] = get_pyfg_data(pyfg_fpath);
    [measurements, var_idx_mapping, dim, true_vals] = get_pyfg_data(pyfg_fpath);
    load_data_time = toc;
    fprintf('Finished loading data in %f seconds.\n', load_data_time);

    % [pose_measurements, range_measurements, var_idx_mapping, dim, true_rots, true_pose_translations, true_landmarks] = get_pyfg_data(pyfg_fpath);
    pose_measurements = measurements.pose_measurements;
    range_measurements = measurements.range_measurements;
    problem.num_pose_priors = length(measurements.pose_priors);
    problem.num_landmark_priors = length(measurements.landmark_priors);

    true_rots = true_vals.true_rots;
    true_pose_translations = true_vals.true_pose_translations;
    true_landmarks = true_vals.true_landmarks;

    %%%%% Useful metadata %%%%%
    problem.dim = dim;
    problem.num_robots = get_num_robots_from_var_idx_map(var_idx_mapping);
    problem.num_poses = var_idx_mapping.pose_var_name_to_cnt.length;
    problem.num_rel_pose_measurements = size(pose_measurements.edges, 1);
    problem.num_landmarks = var_idx_mapping.landmark_name_to_cnt.length;
    problem.num_range_measurements = size(range_measurements.edges, 1);

    % establish variable indices
    problem.all_R_idxs = 1:problem.num_poses*dim;
    last_rot_idx = problem.all_R_idxs(end);
    problem.all_d_idxs = last_rot_idx + (1:problem.num_range_measurements);
    last_d_idx = problem.all_d_idxs(end);
    problem.all_t_idxs = last_d_idx + (1:problem.num_poses);
    last_t_idx = problem.all_t_idxs(end);
    problem.all_l_idxs = last_t_idx + (1:problem.num_landmarks);

    % establish X_gt
    true_dists = get_true_dist_vectors(range_measurements, true_pose_translations, true_landmarks, var_idx_mapping);
    problem.X_gt = [true_rots, true_dists, true_pose_translations, true_landmarks]';

    % establish X_odom
    [odom_rotations, odom_translations] = get_X_odom(pose_measurements, true_vals, var_idx_mapping);
    problem.X_odom = problem.X_gt;
    problem.X_odom(problem.all_R_idxs, :) = odom_rotations';
    problem.X_odom(problem.all_t_idxs, :) = odom_translations';

    %%%%% Pose Graph Component %%%%%
    tic;
    Lrho = construct_connection_Laplacian(pose_measurements, var_idx_mapping);
    assert(size(Lrho, 1) == problem.num_poses * dim);
    assert(size(Lrho, 2) == problem.num_poses * dim);
    con_lap_time = toc;
    fprintf('Finished constructing connection Laplacian in %f seconds.\n', con_lap_time);

    tic;
    Apose = construct_pose_incidence_matrix(pose_measurements, var_idx_mapping);
    assert(size(Apose, 1) == problem.num_rel_pose_measurements);
    assert(size(Apose, 2) == problem.num_poses);
    incidence_time = toc;
    fprintf('Finished constructing incidence matrix in %f seconds.\n', incidence_time);

    tic;
    [T, Omega] = construct_translational_matrices(pose_measurements, var_idx_mapping);
    assert(size(T, 1) == problem.num_rel_pose_measurements);
    assert(size(T, 2) == problem.num_poses * dim);
    assert(size(Omega, 1) == problem.num_rel_pose_measurements);
    assert(size(Omega, 2) == problem.num_rel_pose_measurements);
    translational_time = toc;
    fprintf('Finished constructing translational matrices in %f seconds.\n', translational_time);

    tic;
    V = Apose' * Omega * T;
    V_time = toc;
    problem.V = V;
    fprintf('Finished constructing V matrix in %f seconds.\n', V_time);

    problem.Lrho = Lrho;
    problem.Apose = Apose;
    problem.Omega = Omega;
    problem.T = T;

    %%%%% Range Component %%%%%
    tic;
    [Arange, W, D] = construct_range_graph_matrices(range_measurements, var_idx_mapping);
    assert(size(Arange, 1) == problem.num_range_measurements);
    assert(size(Arange, 2) == problem.num_poses + problem.num_landmarks);
    assert(size(W, 1) == problem.num_range_measurements);
    assert(size(W, 2) == problem.num_range_measurements);
    assert(size(D, 1) == problem.num_range_measurements);
    assert(size(D, 2) == problem.num_range_measurements);
    incidence_time = toc;
    fprintf('Finished constructing range matrices A, D, and W in %f seconds.\n', incidence_time);

    problem.Arange = Arange;
    problem.W = W;
    problem.D = D;

    %%%%% The main Q matrix %%%%%
    Sigma = T' * Omega * T;
    problem.Sigma = Sigma;

    % range graph laplacian
    Lrange = Arange' * W * Arange;
    problem.Lrange = Lrange;

    % translation graph laplacian, embedded in the range graph laplacian space
    Ltau = Apose' * Omega * Apose;
    Ltau = zeroPadSparseSquareMatrix(Ltau, size(Lrange, 1));
    problem.Ltau = Ltau;

    % some zero matrices for padding
    n = problem.num_poses;
    l = problem.num_landmarks;
    r = problem.num_range_measurements;
    ZeroRange = sparse(dim*n, r);
    ZeroL = sparse(dim*n, l);

    Lrot = Lrho + Sigma;
    Ltrans = Ltau + Lrange;
    paddedVT = [V', ZeroL];
    Q = [Lrot,                      ZeroRange,                   paddedVT;
         ZeroRange',                     W * D * D,              D * W * Arange;
         paddedVT',               Arange' * W * D,        Ltrans];
    problem.Q = Q;

    if use_marginalized

        %{
            Qmarg = [Lrho + Sigma, Zero;]    +  [Q11, Q12;]
                    [Zero',        W*D*D]       [Q12', Q22]

            We call the first matrix Qmain, and the second Qtrans.

            Where the second matrix is:
                Qtrans =
                    [V', ZeroL; D*W*Arange] * inv(Ltrans) * [V', ZeroL; D*W*Arange]'

            To efficiently compute inv(Ltrans) * X, we use the Cholesky decomposition
        %}
        tic;
        problem.Qmain = [   Lrho + Sigma,              ZeroRange,      ;
                            ZeroRange',                     W * D * D ];

        % the product operator and laplacian for the marginalized product
        AposePadded = sparse(size(Apose, 1), size(Arange, 2));
        AposePadded(:, 1:size(Apose, 2)) = Apose;
        problem.LeftOperator = [T' * Omega * AposePadded;
                                D * W * Arange];
        problem.Ltrans = Ltrans;

        % the reduced version of the left operator is the same as the full one, but without the last column
        problem.LeftOperatorRed = problem.LeftOperator(:, 1:end-1);
        problem.LtransCholRed = chol(Ltrans(1:end-1, 1:end-1), 'lower');
        marginalized_prod_time = toc;
        fprintf('Finished constructing matrices for marginalized product in %f seconds.\n', marginalized_prod_time);
    end

    if use_marginalized
        problem.block_size = size(problem.Qmain, 1);
    else
        problem.block_size = size(problem.Q, 1);
    end
    problem.stacked_constraints = get_constraints_as_stacked_block_matrix(problem);


end

function paddedMatrix = zeroPadSparseSquareMatrix(originalMatrix, n)
    % Create an nxn sparse matrix filled with zeros
    paddedMatrix = sparse(n, n);

    % Replace the top-left m x m submatrix with the original matrix
    m = size(originalMatrix, 1);
    paddedMatrix(1:m, 1:m) = originalMatrix;
end

function [rotations, translations] = get_X_odom(pose_measurements, true_vals, var_idx_mapping)
    %{
    Returns the odometry matrix X_odom, which is the initialization obtained by
    composing the odometry measurements (as determined by rel-pose measurements between
    subsequent poses)

    Note: this relies on the fact that PYFG files write the odometry measurements
    before any other relative-pose measurements!!
    %}
    dim = length(pose_measurements.t{1});

    true_pose_rots = true_vals.true_rots;
    true_pose_translations = true_vals.true_pose_translations;
    pose_var_names = var_idx_mapping.pose_var_name_to_cnt.keys;

    rotations = [];
    translations = [];
    num_previous_odom_measures = 0;
    num_robots = get_num_robots_from_var_idx_map(var_idx_mapping);
    robot_names = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"];
    for robot_idx = 1:num_robots
        robot_i_name = robot_names(robot_idx);
        robot_i_pose_vars = pose_var_names(startsWith(pose_var_names, robot_i_name));

        % sort the pose vars by their number; format: <Letter><Number>
        % we want A0, A1, A2, ... A10, A11, ... A100, A101, ...
        % so we want to isolate the number and sort by that
        robot_i_pose_idxs = cellfun(@(x) str2double(x(2:end)), robot_i_pose_vars);
        [~, sort_idxs] = sort(robot_i_pose_idxs);
        robot_i_pose_vars = robot_i_pose_vars(sort_idxs);

        start_pose_idx = var_idx_mapping.pose_var_name_to_cnt(cell2mat(robot_i_pose_vars(1)));

        current_T = eye(dim+1);
        start_rot_idxs = (1:dim) + dim * (start_pose_idx-1);
        start_trans_idx = start_pose_idx;
        current_T(1:dim, 1:dim) = true_pose_rots(:, start_rot_idxs);
        current_T(1:dim, end) = true_pose_translations(:, start_trans_idx);

        rotations = [rotations, current_T(1:dim, 1:dim)];
        translations = [translations, current_T(1:dim, end)];

        num_robot_i_poses = length(robot_i_pose_vars);
        num_odom_measures = num_robot_i_poses - 1;
        odom_measure_idxs = num_previous_odom_measures + (1:num_odom_measures);
        assert(length(odom_measure_idxs) == num_odom_measures);
        for odom_idx = odom_measure_idxs
            odom_measure_t = pose_measurements.t{odom_idx};
            odom_measure_R = pose_measurements.R{odom_idx};
            deltaT = [odom_measure_R, odom_measure_t; zeros(1, dim), 1];

            current_T = current_T * deltaT;
            rotations = [rotations, current_T(1:dim, 1:dim)];
            translations = [translations, current_T(1:dim, end)];

        end
        fprintf('Finished robot %s\n', robot_i_name);
        num_previous_odom_measures = num_previous_odom_measures + num_odom_measures;

    end
    disp("done");

end

function true_dist_vectors = get_true_dist_vectors(range_measurements, true_pose_translations, true_landmarks, var_idx_mapping)

    num_pose_vars = size(true_pose_translations, 2);

    % join the translations and landmarks into a single matrix (d x [n + l])
    joined_true_locations = [true_pose_translations, true_landmarks];

    % create a corresponding mapping from vertex name to index in the joined matrix
    all_keys = [var_idx_mapping.pose_var_name_to_cnt.keys, var_idx_mapping.landmark_name_to_cnt.keys];
    offset_landmark_vals = cellfun(@(x) x + num_pose_vars, var_idx_mapping.landmark_name_to_cnt.values, 'UniformOutput', false);
    all_values = [var_idx_mapping.pose_var_name_to_cnt.values, offset_landmark_vals];
    joined_mapping = containers.Map(all_keys, all_values);

    % get the vertex names for the range measurements and find their indices in the joined matrix
    range_measure_out_vertices = range_measurements.edges(:, 1);
    range_measure_in_vertices = range_measurements.edges(:, 2);
    true_out_idx = cell2mat(values(joined_mapping, range_measure_out_vertices));
    true_in_idx = cell2mat(values(joined_mapping, range_measure_in_vertices));

    % get the positions of the vertices in the joined matrix
    true_out_pos = joined_true_locations(:, true_out_idx);
    true_in_pos = joined_true_locations(:, true_in_idx);

    % get the difference between the positions
    true_dist_vectors = true_in_pos - true_out_pos;

    % normalize each column to be unit length
    true_dist_vectors = true_dist_vectors ./ vecnorm(true_dist_vectors);

end

function num_robots = get_num_robots_from_var_idx_map(var_idx_mapping)
    firstChars = cellfun(@(x) x(1), var_idx_mapping.pose_var_name_to_cnt.keys, 'UniformOutput', false);

    % Convert cell array of characters to character array
    charArray = char(firstChars);

    % Get unique characters
    uniqueChars = unique(charArray);

    % Count the number of unique characters
    num_robots = numel(uniqueChars);
end