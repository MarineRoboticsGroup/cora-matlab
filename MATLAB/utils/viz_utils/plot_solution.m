function plot_solution(X, problem_data)

    % make sure that the height of X is dim
    dim = problem_data.dim;
    if size(X, 1) ~= dim
        X = X';
    end
    assert(size(X, 1) == dim);

    X = align_solution_by_first_pose(X, problem_data);

    % set up the figure for plotting
    figure(1)
    clf;
    hold on;

    %%%%%%%% PLOT POSES

    has_priors = (problem_data.num_pose_priors + problem_data.num_beacon_priors) > 0;

    num_robot_poses = problem_data.num_poses-has_priors; %remove the aux pose
    num_robots = problem_data.num_robots;
    poses_per_robot = num_robot_poses / num_robots;
    robot_poses_idx_length = poses_per_robot * (dim+1);

    if num_robots > 0
        assert(mod(num_robot_poses, num_robots) == 0);

        skip_aux_pose_offset = has_priors * (dim+1);

        robot_pose_start_idx = 1 + skip_aux_pose_offset;
        last_pose_idx = num_robot_poses * (dim+1) + skip_aux_pose_offset;

        all_robot_poses = X(:, robot_pose_start_idx:last_pose_idx);
        trans_idxs = dim+1:dim+1:robot_poses_idx_length; % trans idxs for each robot

        % iterate over robots and plot each trajectory with deterministic coloring
        colors = hsv(num_robots);
        for robot_idx = 1:num_robots
            robot_i_poses_start = (robot_idx-1)*robot_poses_idx_length+1;
            robot_i_poses_end = robot_i_poses_start + robot_poses_idx_length-1;
            robot_i_poses = all_robot_poses(:, robot_i_poses_start:robot_i_poses_end);
            assert(length(robot_i_poses) == robot_poses_idx_length);

            trans_i = robot_i_poses(:, trans_idxs);
            assert(length(trans_i)==poses_per_robot);

            plot(trans_i(1, :), trans_i(2, :), 'Color', colors(robot_idx, :), 'LineWidth', 2);
        end
    end

    %%%%%%%% PLOT LANDMARKS

    num_lands = problem_data.num_landmarks;
    if num_lands > 0
        landmarks = X(:, problem_data.all_l_idxs);
    else
        landmarks = [];
    end

    % scatter plot the landmarks
    if ~isempty(landmarks)
        scatter(landmarks(1, :), landmarks(2, :), 20, 'bo')
    end
    
    vmax = max(X');
    vmin = min(X');

    xlim([vmin(1)-3, vmax(1)+3]);
    ylim([vmin(2)-3, vmax(2)+3]);

    % remind to still try plotting ranges between objects, etc
    warning('You should try plotting ranges between objects, etc')

end


