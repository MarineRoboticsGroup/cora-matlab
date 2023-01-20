function file_paths = get_experiment_filepaths_in_subdirs(search_dir)

    % check that the directory exists
    if ~exist(search_dir, 'dir')
        error('Directory %s does not exist', search_dir);
    end

    file_substr = "factor_graph";

    % recursively search the directory for files that begin with file_substr
    file_paths = dir(fullfile(search_dir, '**', file_substr + "*.mat"));

    % throw out any directories that were found
    file_paths = file_paths(~[file_paths.isdir]);

    % don't keep files that end with results.mat or solver_info.mat
    file_paths = file_paths(~endsWith({file_paths.name}, 'results.mat'));
    file_paths = file_paths(~endsWith({file_paths.name}, 'solver_info.mat'));

    % get the full path of each file as a string
    file_paths = string({file_paths.folder}) + "/" + string({file_paths.name});

end