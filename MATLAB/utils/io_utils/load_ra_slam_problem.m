function problem = load_ra_slam_problem(mats_filepath)
    % Q = problem.Q;
    % num_poses = problem.num_poses;
    % num_landmarks = problem.num_landmarks;
    % num_ranges = problem.num_range_measurements;
    % dim = problem.dim;

    problem = struct();

    % load the data from the file
    mats = load(mats_filepath);

    % copy over everything from mats to problem
    fields = fieldnames(mats);
    for i = 1:length(fields)

        % skip the file if it contains "idxs" - we have to handle those a bit
        % differently (see below)
        if ~isempty(strfind(fields{i}, 'idxs'))
            continue;
        end
        data = mats.(fields{i});

        % if the field begins with "num_", then let's cast it to a double
        if strfind(fields{i}, 'num_') == 1
            data = double(data);
        end

        % if the field is "dim" let's cast it to a double
        if strcmp(fields{i}, 'dim')
            data = double(data);
        end

        problem.(fields{i}) = data;
    end


    % indices of variables
    problem.all_R_idxs = mats.rot_idxs+1;
    problem.all_t_idxs = mats.tran_idxs+1;
    problem.all_l_idxs = mats.beacon_idxs+1;
    problem.all_d_idxs = mats.range_idxs+1;

    % metadata on the problem
    problem.num_robots = double(problem.num_robots);
    problem.num_poses = double(problem.num_poses);
    problem.num_landmarks = double(problem.num_landmarks);
    problem.num_range_measurements = double(problem.num_range_measurements);
    problem.dim = double(mats.dim);


end