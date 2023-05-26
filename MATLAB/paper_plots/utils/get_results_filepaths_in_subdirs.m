function file_pairs = get_results_filepaths_in_subdirs(search_dir)

    % check that the directory exists
    if ~exist(search_dir, 'dir')
        error('Directory %s does not exist', search_dir);
    end

    file_substr = "factor_graph";

    % recursively search the directory for files that begin with file_substr
    file_paths = dir(fullfile(search_dir, '**', file_substr + "*.mat"));

    % throw out any directories that were found
    file_paths = file_paths(~[file_paths.isdir]);

    % only keep files that end with results.mat or solver_info.mat
    file_paths = file_paths(endsWith({file_paths.name}, 'results.mat') | endsWith({file_paths.name}, 'solver_info.mat'));

    % check that there are 2n files
    if mod(length(file_paths), 2) ~= 0
        error('There are an odd number of files in %s', search_dir);
    end

    % some useful constants
    num_file_pairs = length(file_paths) / 2;
    file_pairs = cell(num_file_pairs, 2);
    result_idx = 1;
    solver_info_idx = 2;

    % pair up the results and solver_info files
    for i = 1:length(file_paths)
        fpath = file_paths(i).folder + "/" + file_paths(i).name;
        if endsWith(fpath, 'results.mat')
            result_fpath = fpath;
            solver_info_fpath = strrep(fpath, 'results.mat', 'solver_info.mat');
        elseif endsWith(fpath, 'solver_info.mat')
            result_fpath = strrep(fpath, 'solver_info.mat', 'results.mat');
            solver_info_fpath = fpath;
        else
            error('File %s does not end with results.mat or solver_info.mat', fpath);
        end
        file_pairs{i, result_idx} = result_fpath;
        file_pairs{i, solver_info_idx} = solver_info_fpath;
    end

end