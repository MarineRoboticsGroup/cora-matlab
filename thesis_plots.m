
% quickly plot a well-conditioned and pooorly conditioned surface
% for the condition number plot
close all;
figure(1);
hold on;
set(gcf, 'Position', get(0, 'Screensize'));

% Define the x and y ranges
x = linspace(-1, 1, 100);
y = linspace(-1, 1, 100);
[X, Y] = meshgrid(x, y);

% draw a narrow quadratic surface and a wide quadratic surface to show the difference
Z1 = X.^2 + 0.1*Y.^2 + 1;
Z2 = X.^2 + Y.^2 + .3;

% Plot the narrow surface
subplot(1, 2, 1);
surf(X, Y, Z1);
% title('Well-conditioned Surface');
% xlabel('x');
% ylabel('y');
% zlabel('z');
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 20);
set(gca, 'ZScale', 'log');
xticks([]);
yticks([]);
zticks([]);

shading interp;

% Plot the wide surface
subplot(1, 2, 2);
surf(X, Y, Z2);
% title('Poorly-conditioned Surface');
% xlabel('x');
% ylabel('y');
% zlabel('z');
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 20);
% set(gca, 'ZScale', 'log');

% don't draw lines on the surface
shading interp;

% no axis labels or ticks
xticks([]);
yticks([]);
zticks([]);

% background white
set(gcf, 'Color', 'w');



%%
use_marginalized_solver = true;
verbose = true;
fpath = "/home/alan/experimental_data/tiers/ros2/factor_graph.pyfg"; %TIER
% fpath = "/home/alan/experimental_data/highbay_single_drone/factor_graph.pyfg";
% fpath = "/tmp/pyfg/outfinite_nov1/fp_threshold_-99/rx_threshold_-99/lower_dist_1.0/upper_dist_999.0/range_window_0.5/compose_nth_pose_10/disable_outlier_rejection_/drop_all_robot_robot_ranges_/drop_all_robot_beacon_ranges_/loop_closure_--add_loop_closures/dist_offset_0.0/robot_robot_dist_offset_0.0/final.pyfg";
prob = parse_pyfg_text(fpath, use_marginalized_solver, verbose);
n = size(prob.Qmain, 1);
npose = prob.num_poses;
d = prob.dim;
l = prob.num_range_measurements;
t = prob.num_landmarks+prob.num_poses;
trans_idxs = prob.all_t_idxs;
range_idxs = prob.all_d_idxs;
rot_idxs = prob.all_R_idxs;
close all;

% QmainX = problem_data.Qmain * X;
%
% QtransX = problem_data.LeftOperatorRed * ...
%  (problem_data.LtransCholRed' \ ...
%     (problem_data.LtransCholRed \ ...
%         (problem_data.LeftOperatorRed' * X)));
%
% QX = QmainX - QtransX;

% Extract matrices
Q1 = prob.Qmain;
Lred = prob.LeftOperatorRed;
LchholRed = prob.LtransCholRed;

% Print the sizes of the matrices
fprintf("Q1 size: %dx%d\n", slkize(Q1, 1), size(Q1, 2));
fprintf("Lred size: %dx%d\n", size(Lred, 1), size(Lred, 2));
fprintf("LchholRed size: %dx%d\n", size(LchholRed, 1), size(LchholRed, 2));

% Plot the three matrices side by side with matching y-axis limits
figure(1);
hold on;
set(gcf, 'Position', get(0, 'Screensize'));

% Define the y-axis limits based on the maximum height among the matrices
yheight = max([size(Q1, 1), size(Lred, 1), size(LchholRed, 1)]);

% make the dots bigger
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 20);

% Plot Q1
subplot(1, 2, 1);
spy(Q1);
ylim([-0.5, yheight + 0.5]); % Ensure the y-axis range matches the height
title('');
xlabel('');
xticks([])
yticks([])
ylabel('');
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 20);


% Plot Lred
subplot(1, 2, 2);
spy(Lred);
ylim([0.5, yheight + 0.5]); % Ensure the y-axis range matches the height
title('');
xlabel('');
ylabel('');
xticks([])
yticks([])
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 20);


% % Plot LchholRed
% subplot(1, 3, 3)

% Adjust spacing between subplots for better readability
set(gcf, 'Position', get(0, 'Screensize'));

% set the background color to white
set(gcf, 'Color', 'w');

hold off;

figure(2)
spy(LchholRed);
% title('LchholRed');
xlabel('');
% ylabel('');
xticks([])
yticks([])
set(gcf, 'Color', 'w');
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 20);




%%

big_square_size = max(trans_idxs);
small_square_size = max(range_idxs);
trans_line = min(trans_idxs);
range_line = min(range_idxs);

% figure(1)
% spy(prob.Q);
% hold on;
% xline(trans_line, LineWidth=2, Color='black');
% yline(trans_line, LineWidth=2, Color='black');
% xline(range_line, LineWidth=2, Color='black');
% yline(range_line, LineWidth=2, Color='black');
% title(sprintf("Q size: %d", size(prob.Q, 1)));

% % make the figure for the marginalized data matrix
tic;
I = speye(small_square_size);
Qmarg = Qproduct(I, prob);
fprintf("Time to compute Qproduct: %f\n", toc);

% figure(2)
% spy(Qmarg);
% xline(trans_line, LineWidth=2, Color='black');
% yline(trans_line, LineWidth=2, Color='black');
% xline(range_line, LineWidth=2, Color='black');
% yline(range_line, LineWidth=2, Color='black');
% title(sprintf("Qmarg size: %d", size(Qmarg, 1)));

% % what is the relative size of the two matrices?
% rel_dim = size(Qmarg, 1) / size(prob.Q, 1);
% rel_nnzs = nnz(Qmarg) / nnz(prob.Q);
fprintf("Relative size of Qmarg to Q: %f\n", rel_dim);
fprintf("Relative nnzs of Qmarg to Q: %f\n", rel_nnzs);


% can we plot the spy of the two matrices side by side?
figure(3)
hold on;
set(gcf, 'Position', get(0, 'Screensize'));

% Plot the larger matrix
subplot(1, 2, 1);
spy(prob.Q);
axis([0, big_square_size, 0, big_square_size]); % Set axis limits to match the larger matrix
xline(trans_line, LineWidth=2, Color='black');
yline(trans_line, LineWidth=2, Color='black');
xline(range_line, LineWidth=2, Color='black');
yline(range_line, LineWidth=2, Color='black');
% title should be 'latex(Q) (Original Matrix)'
% title_str = sprintf('$$Q$$ (Original Matrix)');
% title(title_str, 'Interpreter', 'latex');

% Plot the smaller matrix
subplot(1, 2, 2);
spy(Qmarg);
axis([0, big_square_size, 0, big_square_size]); % Match axis limits of the larger matrix for proportion
xline(trans_line, LineWidth=2, Color='black');
yline(trans_line, LineWidth=2, Color='black');
xline(range_line, LineWidth=2, Color='black');
yline(range_line, LineWidth=2, Color='black');
% title_str = sprintf('$$Q_{marg}$$ (Marginalized Matrix)');
% title(title_str, 'Interpreter', 'latex');

% Add a main title with relative dimensions and non-zero entries
% sgtitle(sprintf("Relative size: %d%%, Relative nnzs: %d%%", round(rel_dim*100), round(rel_nnzs*100)));

% set the background color to white
set(gcf, 'Color', 'w');


