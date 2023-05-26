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
    gt_color = [0.5, 0.5, 0.5];

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

    % set linewidth and marker size
    set(fig, 'DefaultLineLineWidth', 2);
    improvePlot();

    % set axis limits
    x_spread = abs(max_x - min_x);
    y_spread = abs(max_y - min_y);
    x_buff = x_spread/5;
    y_buff = y_spread/10;
    xlim([min_x - x_buff, max_x + x_buff]);
    ylim([min_y - y_buff, max_y + y_buff]);

    % make the axis square
    axis square

    % Xvals_rounded is a 3D array of size (num_iterates, num_rows, num_cols)
    % we want to animate the 2D slices of the 3D array by slicing along the
    % first dimension
    Xvals = squeeze(Xvals_rounded(1, :, :));
    num_robots = problem_data.num_robots;

    % make all of the robot plot objects
    fprintf('Generating all of the plots for the robot trajectories\n')
    robot_plots = cell(num_robots, 1);
    robot_t_idxs = cell(num_robots, 1);
    % get a different color for each robot
    colors = lines(num_robots);
    for robot_idx = 1:num_robots
        robot_t_idxs{robot_idx} = get_robot_t_idxs(problem_data, robot_idx);
        if show_gt
            % plot gt trajectory as black dashed line with line width 1
            plot(gt_vals(1, robot_t_idxs{robot_idx}), gt_vals(2, robot_t_idxs{robot_idx}), '-', 'Color', gt_color, 'LineWidth', 1);
        end
        % robot_plots{robot_idx} = plot(Xvals(1, robot_t_idxs{robot_idx}), Xvals(2, robot_t_idxs{robot_idx}), 'Color', "green");
        robot_plots{robot_idx} = plot(Xvals(1, robot_t_idxs{robot_idx}), Xvals(2, robot_t_idxs{robot_idx}), 'Color', colors(robot_idx, :));
    end

    % make all of the landmark plot objects
    fprintf('Generating all of the plots for the landmarks\n')
    if show_gt
        % plot gt_vals as grey x's
        scatter(gt_vals(1, l_idxs), gt_vals(2, l_idxs), 30, 'x', 'MarkerEdgeColor', gt_color, 'LineWidth', 2);
    end
    l = scatter(Xvals(1,problem_data.all_l_idxs), Xvals(2, problem_data.all_l_idxs), 30, 'bo', "LineWidth", 2);

    % set a legend if we are showing the ground truth
    if show_gt
        gt_name = "ground truth";
        est_name = "CORA";
        legend(gt_name, est_name, gt_name, est_name, "Location","northwest");
    end

    % run over the iterates and show all of the plot
    fprintf('Performing plotting and saving frames to make .gif\n')
    gif_fpath = strrep(data_path, ".mat", "_projected_iterates.gif");
    video_fpath = strrep(data_path, "factor_graph.mat", "cora_animation.avi");
    make_new_gif = true;
    save_images = false;

    if make_new_gif
        im = cell(num_iterates, 1);
        vidWriter = VideoWriter(video_fpath);
        open(vidWriter);
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

        if save_images
            img_dir = strrep(data_path, ".mat", "_projected_iterates");
            if ~exist(img_dir, 'dir')
                mkdir(img_dir);
            end

            % make img path, with idx padded to 4 digits
            img_path = fullfile(img_dir, sprintf('img_%04d.png', idx));

            % write the image to file with a high resolution
            resolution = 600;
            exportgraphics(fig, img_path, "Resolution", resolution);

            % save the image with a high dpi
            % print(fig, img_path, '-dpng', '-r300');
            % print(fig, img_path, '-dpng', '-r600');

            fprintf('Saved image to %s\n', img_path);
        end

        if make_new_gif
            im{idx} = frame2im(frame);
            writeVideo(vidWriter, im{idx});
        end
    end

    if  make_new_gif
        delay_between_frames_sec = 0.01;
        fprintf('Saving gif to %s\n', gif_fpath);
        for idx = 1:num_iterates
            [A,map] = rgb2ind(im{idx},128);
            if idx == 1
                imwrite(A,map,gif_fpath,"gif","LoopCount",Inf,"DelayTime",delay_between_frames_sec);
            else
                imwrite(A,map,gif_fpath,"gif","WriteMode","append","DelayTime",delay_between_frames_sec);
            end
        end
        fprintf('Saved gif to %s\n', gif_fpath);
        fprintf('Saved movie to %s\n', video_fpath);
        close(vidWriter);
    end
end

function [] = improvePlot()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot parameters
    % MATLAB treats mac and PC displays differently which can create
    % weird looking graphs. Here we handle system differences

    if ismac
        plot_width_in_px = 800;
        plot_height_in_px = 800;
        marker_size=15;
        marker_line_width=2.5;
        box_thickness = 3;
        axis_tick_font_size = 24;
        axis_label_font_size = 24;
        legend_font_size = 20;
        error_bar_cap_size = 15;
    else % (ispc || isunix)
        plot_width_in_px = 600;
        plot_height_in_px = 600;
        marker_size=10;
        marker_line_width=2.0;
        box_thickness = 2;
        axis_tick_font_size = 18;
        axis_label_font_size = 18;
        legend_font_size = 16;
        error_bar_cap_size = 10;
    end

    marker_outline = 'matching'; % could be 'black' or 'matching'

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Use h as handle for current figure
    hFig = gcf;
    % Change figure background colour to white
    set(hFig, 'Color', 'white');

    % Make the figure bigger
    set(hFig, 'rend', 'painters', 'Units', 'pixels', 'pos', ...
        [100 100 plot_width_in_px plot_height_in_px]);

    % Grab the axes handle(s)
    axis_handles=findobj(hFig,'type','axe');

    % Iterate over all axes handle(s), this is useful if there are subplots
    for i = 1:length(axis_handles)
        ax = axis_handles(i);

        % Change default font size (tick labels, legend, etc.)
        set(ax, 'FontSize', axis_tick_font_size, 'FontName', 'Arial', 'LineWidth', box_thickness);

        set(ax, 'Box', 'on');

        % Change font size for axis text labels
        set(get(ax, 'XLabel'),'FontSize', axis_label_font_size, 'FontWeight', 'Bold');
        set(get(ax, 'YLabel'),'FontSize', axis_label_font_size, 'FontWeight', 'Bold');

        try % try statement to avoid error with categorical axes
        ax.XRuler.Exponent = 0; % Remove exponential notation from the X axis
        ax.YRuler.Exponent = 0; % Remove exponential notation from the Y axis
        catch
        end

    end

    % Find all the lines, and markers
    LineH = findobj(hFig, 'type', 'line', '-or', 'type', 'errorbar');

    if(~isempty(LineH))
        for i=1:length(LineH) % Iterate over all lines in the plot
            % Decide what color for the marker edges
            this_line_color = get(LineH(i),'color');
            if strcmp(marker_outline, 'black')
                marker_outline_color = 'black';
            elseif strcmp(marker_outline, 'matching')
                marker_outline_color = this_line_color;
            else
                marker_outline_color = 'black';
            end

            % If the LineWidth has not been customized, then change it
            if (get(LineH(i), 'LineWidth') <= 1.0)
                set(LineH(i), 'LineWidth', marker_line_width)
            end
            % Change lines and markers if they exist on the plot
            set(LineH(i),   'MarkerSize', marker_size, ...
                'MarkerEdgeColor', marker_outline_color, ...
                'MarkerFaceColor', this_line_color);
        end
    end

    % Find and change the error bars
    LineH = findobj(hFig, 'type', 'errorbar');
    if(~isempty(LineH))
        for i=1:length(LineH) % Iterate over all lines in the plot
            LineH(i).CapSize=error_bar_cap_size;
%             LineH(i).Color = [0 0 0]; % Set all error bars to black

        end
    end

    % Find the legend, and if there is one, change it
    h = get(hFig,'children');
    for k = 1:length(h)
        if strcmpi(get(h(k),'Tag'),'legend')
            set(h(k), 'FontSize', legend_font_size, 'location', 'northwest');
            break;
        end
    end

end