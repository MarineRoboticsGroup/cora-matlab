clc; clear all; close all;

% whether loop closures or not
no_loop_base_data_dir = "~/data/manhattan/cert/no_loop_closures/";
yes_loop_base_data_dir = "~/data/manhattan/cert/100loop_closures/";
data_dir_options = [no_loop_base_data_dir, yes_loop_base_data_dir];
num_data_dirs = length(data_dir_options);

% different param sweeps
sweep_options = ["sweep_num_poses", "sweep_num_ranges",  "sweep_num_robots",  "sweep_range_cov"];
num_options = length(sweep_options);

% where to save the figure to
repo_fig_base_dir = "~/rss23-ra-slam-certification/figures/data_analysis/";
save_data_dir = repo_fig_base_dir;


fig = clf(figure(3));
fig.WindowState = 'maximized';

for data_dir_idx = 1:num_data_dirs
    base_data_dir = data_dir_options(data_dir_idx);
    for options_idx = 1:num_options
        sweep_opt = sweep_options(options_idx);
        sweep_data_dir = base_data_dir + sweep_opt + "/";
        results_fpaths_pairs = get_results_filepaths_in_subdirs(sweep_data_dir);

        if sweep_opt == "sweep_num_poses"
            param = get_num_poses_from_fpath_names(results_fpaths_pairs);
            param_name = "# poses";
        elseif sweep_opt == "sweep_num_ranges"
            param = get_num_ranges_from_fpath_names(results_fpaths_pairs);
            param_name = "# range measurements";
        elseif sweep_opt == "sweep_num_robots"
            param = get_num_robots_from_fpath_names(results_fpaths_pairs);
            disp(param);
            param_name = "# robots";
        elseif sweep_opt == "sweep_range_cov"
            param = get_range_covs_from_fpath_names(results_fpaths_pairs);
            % this param name should be sigma^2 in the fancy latex notation
            param_name = "\sigma_{ij}^2";
        end

        % get the suboptimality bounds
        [abs_gaps, rel_gaps] = get_suboptimality_gaps_from_results_fpaths(results_fpaths_pairs);

        % convert rel_gaps to percent
        rel_gaps = rel_gaps * 100;

        % combine all of the information into a single table
        results_table = table(param, abs_gaps, rel_gaps);
        results_table = sortrows(results_table, 'param');

        % plot the relative suboptimality gap
        plot_idx = (data_dir_idx - 1) * num_options + options_idx;
        subplot(num_data_dirs,num_options,plot_idx);

        % plot the param and rel gap with a line and closed circles, linewidth 2
        plot(results_table.param, results_table.rel_gaps, '-o', 'LineWidth', 2);

        % xlabel with fontsize 14
        xlabel(param_name, 'FontSize', 14);
        ylabel("relative suboptimality gap (%)", 'FontSize', 14);

        if data_dir_idx == 1
            ylim([0, 30]);
        else
            ylim([0, 1]);
        end
        xlim([0, max(results_table.param) + 1]);


        % save the plot

    end
end

subplotLayout = [num_data_dirs, num_options];
ax = gobjects(fliplr(subplotLayout));
for i = 1:prod(subplotLayout)
    ax(i) = subplot(subplotLayout(1),subplotLayout(2),i);
end
set(ax,'Units','Normalize') % this is typically the default
ax = ax'; % now axis handles are same shape as subplot layout
% Reduce vertical space of axes just a bit to make room for "titles"
axPos = cell2mat(get(ax, 'Position'));
axPos(:,4) = axPos(:,4).*.96; % reduces vert space to 96% of height
set(ax, {'Position'}, mat2cell(axPos, ones(numel(ax),1), 4))
% Get upper position of each row of axes, normalized coordinates
axPos = cell2mat(get(ax(:,1), 'Position'));
axUpperPos = sum(axPos(:,[2,4]),2);  %upper pos.
% Get center position for 1st row (assumes all rows have same center)
axPos = cell2mat(get(ax(1,[1,end]),'Position'));
axCenterPos = mean([axPos(1,1), sum(axPos(2,[1,3]))]);
% list all titles in order (top to bottom)
titles = {'0 interrobot pose measurements', '100 interrobot pose measurements'};
% Set annotation for each row of subplots
titleHandles = gobjects(numel(titles),1);
for i = 1:numel(titles)
    titleHandles = annotation('textbox','String',titles{i}, ...
        'Position', [axCenterPos, axUpperPos(i), 0, 0], ...
        'HorizontalAlignment', 'center','VerticalAlignment','bottom',...
        'LineStyle','none','FitBoxToText','on', ...
        'FontWeight',ax(1).Title.FontWeight, ... % matches title property
        'FontSize', 18, ...    % matches title property
        'FontName', ax(1).Title.FontName, ...    % matches title property
        'Color', ax(1).Title.Color);             % matches title property
end

plot_fpath = save_data_dir + "suboptimality_gaps_vs_params.png";
saveas(fig, plot_fpath);
fprintf("Saved plot to %s\n", plot_fpath);
