% function [pose_measurements, range_measurements, var_idx_mapping, d, true_rots, true_pose_translations, true_landmarks] = get_pyfg_data(pyfg_fpath)
function [measurements, var_idx_mapping, d, true_vals] = get_pyfg_data(pyfg_fpath)
    %
    % a 2D or 3D pose graph SLAM problem, and returns a MATLAB struct
    % 'pose_measurements' containing the description of the problem in the format
    % required by SE-Sync.  More precisely, 'pose_measurements' consists of:
    %
    % edges:  An (mx2)-dimensional matrix encoding the edges in the measurement
    %      network; edges(k, :) = [i,j] means that the kth measurement is of the
    %      relative transform x_i^{-1} x_j.  NB:  This indexing scheme requires
    %      that the states x_i are numbered sequentially as x_1, ... x_n.
    % R:  An m-dimensional cell array whose kth element is the rotational part
    %      of the kth measurement
    % t:  An m-dimensional cell array whose kth element is the translational
    %      part of the kth measurement
    % kappa:  An m-dimensional cell array whose kth element gives the
    %      precision of the rotational part of the kth measurement.
    % tau:  An m-dimensional cell array whose kth element gives the precision
    %      of the translational part of the kth measurement.
    %
    %


    fid = fopen(pyfg_fpath, 'r');

    % variables
    pose_var_names = strings(0);
    landmark_names = strings(0);

    % pose priors
    pose_prior_names = strings(0);
    pose_prior_means = [];
    pose_prior_precisions = [];

    % landmark priors
    landmark_prior_names = strings(0);
    landmark_prior_means = [];
    landmark_prior_precisions = [];

    % measurements
    pose_edge_id = 0;
    range_edge_id = 0;

    % ground truth
    true_rots = [];
    true_pose_translations = [];
    true_landmarks = [];

    read_line = fgets(fid);  % Read the next line from the file

    while ischar(read_line)  % As long as this line is a valid character string
        token = strtok(read_line);

        if (strcmp(token, 'VERTEX_SE3:QUAT'))
            assert(~exist('d', 'var') || d == 3, "3D Pose in 2D problem")
            d = 3;
            % 3D POSE
            % VERTEX_SE3:QUAT <timestamp> <VertexName> <x> <y> <z> <qx> <qy> <qz> <qw>
            % VERTEX_SE3:QUAT 1691505493.647804737 A1 5.000000 0.000000 -10.000000 0.0000000 0.0000000 0.0000000 1.0000000
            C = textscan(read_line, '%s %f %s %f %f %f %f %f %f %f');
            [~, ~, var_name, px, py, pz, qx, qy, qz, qw] = C{:};
            pose_var_names = [pose_var_names; var_name];

            true_rots = [true_rots, quat2rotm([qw, qx, qy, qz])];
            true_pose_translations = [true_pose_translations, [px; py; pz]];

        elseif(strcmp(token, 'VERTEX_SE2'))
            assert(~exist('d', 'var') || d == 2, "2D Pose in 3D problem")
            d = 2;
            % 2D POSE
            % VERTEX_SE2 <timestamp> <VertexName> <x> <y> <theta>
            % VERTEX_SE2 3152.0 A0 -34.2086489999201 45.300763999112 1.1205036535897932
            C = textscan(read_line, '%s %f %s %f %f %f');
            [~, ~, var_name, px, py, ptheta] = C{:};
            pose_var_names = [pose_var_names; var_name];

            % Store the true pose
            true_rots = [true_rots, get_rot_mat_from_theta(ptheta)];
            true_pose_translations = [true_pose_translations, [px; py]];

        elseif(strcmp(token, 'VERTEX_XY'))
            assert(~exist('d', 'var') || d == 2, "2D landmark in 3D problem")
            d = 2;
            % 2D LANDMARK
            % VERTEX_XY <VertexName> <x> <y>
            % VERTEX_XY L3 1.70946300006472 -5.81220300029963
            C = textscan(read_line, '%s %s %f %f');
            [~, var_name, lx, ly] = C{:};
            landmark_names = [landmark_names; var_name];
            true_landmarks = [true_landmarks, [lx; ly]];

        elseif(strcmp(token, 'VERTEX_XYZ'))
            % d should either be empty or 3
            assert(~exist('d', 'var') || d == 3, "3D landmark in 2D problem")
            d = 3;
            % 3D LANDMARK
            % VERTEX_XYZ <VertexName> <x> <y> <z>
            % VERTEX_XYZ L3 1.70946300006472 -5.81220300029963 0.0
            C = textscan(read_line, '%s %s %f %f %f');
            [~, var_name, lx, ly, lz] = C{:};
            landmark_names = [landmark_names; var_name];
            true_landmarks = [true_landmarks, [lx; ly; lz]];
        elseif(strcmp(token, 'VERTEX_XY:PRIOR'))
            assert(~exist('d', 'var') || d == 2, "2D landmark prior in 3D problem")
            d = 2;
            % VERTEX_XY:PRIOR 1665936917 L0 20.742063186872322 -25.000484838943855 0.20000000298023224 0.0 0.20000000298023224
            % "{pose_prior_type} {pose_prior.timestamp:.{time_fprec}f} {pose_prior.name} {measurement_values} {measurement_noise}"
            C = textscan(read_line, '%s %f %s %f %f %f %f %f');
            [~, ~, var_name, lx, ly, C11, C12, C22] = C{:};

            % right now assume that the prior is spherical
            assert(C11 == C22 || C12 == 0, "Covariance matrix is not spherical")

            % make sure we haven't already seen this landmark
            assert(~ismember(var_name, landmark_prior_names), "Duplicate landmark prior")

            % Store the landmark prior
            landmark_prior_names = [landmark_prior_names; var_name];
            landmark_prior_means = [landmark_prior_means, [lx; ly]];
            landmark_prior_precisions = [landmark_prior_precisions, 1/C11];

        elseif(strcmp(token, 'EDGE_SE3:QUAT'))
            % 3D Measurement
            assert (d == 3, "3D observation in 2D problem")

            pose_edge_id = pose_edge_id + 1;  % Increment the count for the number of edges

            % The pyfg format specifies a 3D relative pose measurement in the
            % following form:

            % EDGE_SE3:QUAT timestamp id1 id2 dx dy dz dqx dqy dqz dqw
            % C11 C12 C13 C14 C15 C16
            %     C22 C23 C24 C25 C26
            %         C33 C34 C35 C36
            %             C44 C45 C46
            %                 C55 C56
            %                     C66

            % strread(read_line, '%s %s %s %f %f %f %f %f %f %f    %f %f %f %f %f %f    %f %f %f %f %f   %f %f %f %f   %f %f %f    %f %f    %f');
            C = textscan(read_line, '%s %f %s %s %f %f %f %f %f %f %f    %f %f %f %f %f %f    %f %f %f %f %f   %f %f %f %f   %f %f %f    %f %f    %f');
            [~, ~, id1, id2, dx, dy, dz, dqx, dqy, dqz, dqw, ...
                C11, C12, C13, C14, C15, C16, ...
                C22, C23, C24, C25, C26, ...
                C33, C34, C35, C36, ...
                C44, C45, C46, ...
                C55, C56, ...
                C66] = C{:};

            % Store the connectivity of this edge
            pose_edges(pose_edge_id, :) = [id1 , id2 ];

            % Store the translational measurement
            t{pose_edge_id} = [dx, dy, dz]';

            % Reconstruct quaternion for relative measurement
            q = [dqw, dqx, dqy, dqz]';
            q = q / norm(q);  % Make sure that this is properly normalized

            % Compute and store corresponding rotation matrix
            R{pose_edge_id} = quat2rotm(q');

            % Reconstruct the information matrix
            measurement_covar = ...
                [C11, C12, C13, C14, C15, C16,
                C12, C22, C23, C24, C25, C26,
                C13, C23, C33, C34, C35, C36,
                C14, C24, C34, C44, C45, C46,
                C15, C25, C35, C45, C55, C56,
                C16, C26, C36, C46, C56, C66];

            % Compute and store the optimal (information-divergence-minimizing)
            % value of the parameter tau

            tau{pose_edge_id} = 3 / trace(measurement_covar(1:3, 1:3) );

            % Extract and store the optimal (information-divergence-minimizing)
            % value of the parameter kappa
            kappa{pose_edge_id} = 3 / (2 *trace(measurement_covar(4:6, 4:6)));

        elseif(strcmp(token, 'EDGE_SE2'))
            % 2D Measurement
            assert (d == 2, "2D observation in 3D problem")

            pose_edge_id = pose_edge_id + 1;

            % The g2o format specifies a 3D relative pose measurement in the
            % following form:

            % EDGE_SE2 timestamp id1 id2 dx dy dtheta, C11, C12, C13, C22, C23, C33
            % EDGE_SE2 3511.29851794243 A3587 A3588 0.332397242830102 0.0 -0.0319438761202058 0.010000000000000002 0.0 -0.0 0.010000000000000002 -0.0 0.0001
            C = textscan(read_line, '%s %f %s %s %f %f %f %f %f %f %f %f %f');
            [~, ~, id1, id2, dx, dy, dth, C11, C12, C13, C22, C23, C33] = C{:};

            % Store the connectivity of this edge
            pose_edges(pose_edge_id, :) = [id1 , id2 ];

            % Store the translational measurement
            t{pose_edge_id} = [dx, dy]';

            % Reconstruct and store the rotational measurement
            R{pose_edge_id} = [cos(dth), -sin(dth);
                        sin(dth), cos(dth)];

            % Reconstruct the information matrix
            measurement_covar = ...
                [C11, C12, C13;
                C12, C22, C23;
                C13, C23, C33];

            % Extract and store an outer approximation for the translational
            % measurement precision
            tau{pose_edge_id} = 2 / trace( measurement_covar(1:2, 1:2));

            % Extract and store an outer approximation for the rotational
            % measurement precision
            kappa{pose_edge_id} = 1/C33;

        elseif(strcmp(token, 'EDGE_RANGE'))
            % range measurement
            range_edge_id = range_edge_id + 1;

            % EDGE_RANGE 3542.90110708028 A3904 L3 26.69216622156613 0.2655983804002221
            C = textscan(read_line, '%s %f %s %s %f %f');
            [~, ~, id1, id2, range, variance] = C{:};

            % Store the connectivity of this edge
            range_edges(range_edge_id, :) = [id1 , id2 ];

            % Store the range measurement
            ranges{range_edge_id} = range;

            % Store the range measurement precision
            variance = variance * 5;
            range_precision{range_edge_id} = 1 / variance;


        else
            error('Unrecognized token in file: %s', token);

        end

        read_line = fgets(fid);
    end


    % Construct and return pose_measurements struct
    pose_measurements.edges = pose_edges;
    pose_measurements.R = R;
    pose_measurements.t = t;
    pose_measurements.kappa = kappa;
    pose_measurements.tau = tau;

    % Construct and return range_measurements struct
    range_measurements.edges = range_edges;
    range_measurements.range = ranges;
    range_measurements.precision = range_precision;

    % Construct and return landmark_priors struct
    landmark_priors.names = landmark_prior_names;
    landmark_priors.mean = landmark_prior_means;
    landmark_priors.precision = landmark_prior_precisions;

    % Construct and return pose_priors struct
    pose_priors.names = pose_prior_names;
    pose_priors.mean = pose_prior_means;
    pose_priors.precision = pose_prior_precisions;

    % Construct and return var_idx_mapping struct
    num_pose_vars = length(pose_var_names);
    num_landmarks = length(landmark_names);

    % fill measurements struct
    measurements.pose_measurements = pose_measurements;
    measurements.range_measurements = range_measurements;

    if ~exist('pose_priors', 'var')
        measurements.pose_priors = [];
    else
        measurements.pose_priors = pose_priors;
    end

    if ~exist('landmark_priors', 'var')
        measurements.landmark_priors = [];
    else
        measurements.landmark_priors = landmark_priors;
    end

    % fill true_vals struct
    true_vals.true_rots = true_rots;
    true_vals.true_pose_translations = true_pose_translations;
    true_vals.true_landmarks = true_landmarks;

    % make sure none of the names are repeated between the two sets
    assert(isempty(intersect(pose_var_names, landmark_names)), "Variable names are not unique")


    % set up the mapping from variable names to indices
    if isempty(pose_var_names)
        var_idx_mapping.pose_var_name_to_cnt = containers.Map();
    else
        var_idx_mapping.pose_var_name_to_cnt = containers.Map(pose_var_names, 1:num_pose_vars);
    end

    if isempty(landmark_names)
        var_idx_mapping.landmark_name_to_cnt = containers.Map();
    else
        var_idx_mapping.landmark_name_to_cnt = containers.Map(landmark_names, 1:num_landmarks);
    end



end

function rot_mat = get_rot_mat_from_theta(theta)
    rot_mat = [cos(theta), -sin(theta);
                sin(theta), cos(theta)];
end