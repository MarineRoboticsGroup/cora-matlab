function animate_lifted_solver_trajectory(data_path, show_gt)
    res_path = strrep(data_path, ".mat", "_results.mat");
    cora_iterates_info_path = strrep(data_path, '.mat', '_cora_iterates_info.mat');

    problem_data = load_ra_slam_problem(data_path);
    load(res_path);
    load(cora_iterates_info_path);

    num_iterates = length(cora_iterates_info);
    [num_cols, num_rows] = size(problem_data.X_odom);
    Xvals_rounded = zeros(num_iterates, num_rows, num_cols);

    gt_vals = align_solution_by_first_pose(problem_data.X_gt', problem_data);
    % the Plaza2 data needs to be flipped 180 degrees to visualize properly
    % rot_angle_rad = deg2rad(180);
    % full_rot = [cos(rot_angle_rad), -sin(rot_angle_rad);
    %             sin(rot_angle_rad), cos(rot_angle_rad)];
    % gt_vals = full_rot* gt_vals;

    for idx = 1:num_iterates
        % for each iterate, we want to round the solution to SE(d) and
        % align it such that the first pose is at the origin
        verbosity = 0;
        Xvals_rounded(idx, :, :) = align_solution_by_first_pose(...
            round_solution(cora_iterates_info(idx).Xvals',problem_data, verbosity), problem_data...
        );
    end

    % get the max and min x and y values for the plot
    t_idxs = problem_data.all_t_idxs;
    l_idxs = problem_data.all_l_idxs;
    all_xval_translations = Xvals_rounded(:, :, t_idxs);
    all_xval_landmarks = Xvals_rounded(:, :, l_idxs);
    max_xvals_translations = max(all_xval_translations, [], [1,3]);
    min_xvals_translations = min(all_xval_translations, [], [1,3]);

    if isempty(all_xval_landmarks)
        max_x = max_xvals_translations(1);
        max_y = max_xvals_translations(2);
        min_x = min_xvals_translations(1);
        min_y = min_xvals_translations(2);
    else
        max_xvals_landmarks = max(all_xval_landmarks);
        max_x = max(max_xvals_landmarks(1), max_xvals_translations(1));
        max_y = max(max_xvals_landmarks(2), max_xvals_translations(2));
        min_xvals_landmarks = min(all_xval_landmarks);
        min_x = min(min_xvals_landmarks(1), min_xvals_translations(1));
        min_y = min(min_xvals_landmarks(2), min_xvals_translations(2));
    end

    fig = clf(figure(2));
    hold on;
    set(fig, 'DefaultLineLineWidth', 2);

    % set axis limits
    x_spread = abs(max_x - min_x);
    y_spread = abs(max_y - min_y);
    xlim([min_x - x_spread/10, max_x + x_spread/10]);
    ylim([min_y - y_spread/10, max_y + y_spread/10]);

    % Xvals_rounded is a 3D array of size (num_iterates, num_rows, num_cols)
    % we want to animate the 2D slices of the 3D array by slicing along the
    % first dimension
    Xvals = squeeze(Xvals_rounded(1, :, :));
    num_robots = problem_data.num_robots;

    % make all of the robot plot objects
    fprintf('Generating all of the plots for the robot trajectories\n')
    robot_plots = cell(num_robots, 1);
    robot_t_idxs = cell(num_robots, 1);
    for robot_idx = 1:num_robots
        robot_t_idxs{robot_idx} = get_robot_t_idxs(problem_data, robot_idx);
        robot_plots{robot_idx} = plot(Xvals(1, robot_t_idxs{robot_idx}), Xvals(2, robot_t_idxs{robot_idx}));
        if show_gt
            plot(gt_vals(1, robot_t_idxs{robot_idx}), gt_vals(2, robot_t_idxs{robot_idx}));
        end
    end

    % make all of the landmark plot objects
    fprintf('Generating all of the plots for the landmarks\n')
    l = scatter(Xvals(1,problem_data.all_l_idxs), Xvals(2, problem_data.all_l_idxs), 20, 'bo');
    if show_gt
        scatter(gt_vals(1, l_idxs), gt_vals(2, l_idxs), 20, 'rx');
    end

    % set a legend if we are showing the ground truth
    if show_gt
        legend("estimate", "ground truth", "estimate", "ground truth");
    end

    % run over the iterates and show all of the plot
    fprintf('Performing plotting and saving frames to make .gif\n')
    gif_fpath = strrep(data_path, ".mat", "_projected_iterates.gif");
    gif_already_exists = exist(gif_fpath, 'file');
    make_new_gif = ~gif_already_exists;
    make_new_gif = true;
    % if the gif exists ask the user if they want to overwrite it
%     if gif_already_exists
%         overwrite = input(sprintf("%s already exists. Overwrite? (y/n): ", gif_fpath), 's');
%         if overwrite == 'y'
%             make_new_gif = true;
%         end
%     end

    if make_new_gif
        im = cell(num_iterates, 1);
    end
    for idx = 1:num_iterates

        Xvals = squeeze(Xvals_rounded(idx, :, :));

        for robot_idx = 1:num_robots
            robot_plots{robot_idx}.XData = Xvals(1, robot_t_idxs{robot_idx});
            robot_plots{robot_idx}.YData = Xvals(2, robot_t_idxs{robot_idx});
        end
        l.XData = Xvals(1,l_idxs);
        l.YData = Xvals(2,l_idxs);
        frame = getframe(fig);
        if make_new_gif
            im{idx} = frame2im(frame);
        end
    end

    if  make_new_gif
        delay_between_frames_sec = 0.01;
        fprintf('Saving gif to %s\n', gif_fpath);
        for idx = 1:num_iterates
            [A,map] = rgb2ind(im{idx},256);
            if idx == 1
                imwrite(A,map,gif_fpath,"gif","LoopCount",Inf,"DelayTime",delay_between_frames_sec);
            else
                imwrite(A,map,gif_fpath,"gif","WriteMode","append","DelayTime",delay_between_frames_sec);
            end
        end
        fprintf('Saved gif to %s\n', gif_fpath);
    end
end