function plot_solution(X, problem_data, show_gt, flip_gt)

% make sure that the height of X is dim
dim = problem_data.dim;
num_robots = problem_data.num_robots;
if size(X, 1) ~= dim
    X = X';
end
assert(size(X, 1) == dim);

if show_gt
    gt_vals = align_solution_by_first_pose(problem_data.X_gt', problem_data);
    % the Plaza2 data needs to be rotated 180 degrees to visualize properly
    if ~exist('flip_gt', 'var')
        flip_gt = false;
    end
    if flip_gt
        rot_angle_rad = deg2rad(180);
        full_rot = [cos(rot_angle_rad), -sin(rot_angle_rad);
                    sin(rot_angle_rad), cos(rot_angle_rad)];
        gt_vals = full_rot* gt_vals;
    end
end
Xvals = align_solution_by_first_pose(X, problem_data);

% set up the figure for plotting
figure(1)
clf;
hold on;

%%%%%%%% PLOT POSES

% iterate over robots and plot each trajectory with deterministic coloring
colors = hsv(num_robots*2);
for robot_idx = 1:num_robots
    robot_t_idxs = get_robot_t_idxs(problem_data, robot_idx);
    plot(Xvals(1, robot_t_idxs), Xvals(2, robot_t_idxs), 'Color', colors(robot_idx, :), 'LineWidth', 2);
    if show_gt
        plot(gt_vals(1, robot_t_idxs), gt_vals(2, robot_t_idxs),...
            "Color", colors(robot_idx+1, :), 'LineWidth', 1);
    end
end

%%%%%%%% PLOT LANDMARKS

num_lands = problem_data.num_landmarks;
if num_lands > 0
    landmarks = Xvals(:, problem_data.all_l_idxs);
    scatter(landmarks(1, :), landmarks(2, :), 20, 'bo')
    if show_gt
        scatter(gt_vals(1, problem_data.all_l_idxs), gt_vals(2, problem_data.all_l_idxs), 20, 'ro')
    end
end

vmax = max(X');
vmin = min(X');

xlim([vmin(1)-3, vmax(1)+3]);
ylim([vmin(2)-3, vmax(2)+3]);

% wait for user to close the figure
waitfor(figure(1));

% remind to still try plotting ranges between objects, etc
%     warning('You should try plotting ranges between objects, etc')

end


