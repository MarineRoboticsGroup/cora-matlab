function R_idxs = get_robot_R_idxs(problem_data, robot_idx)

    % has_priors = (problem_data.num_pose_priors + problem_data.num_landmark_priors) > 0;
    warning("Assuming no priors - get_robot_R_idxs")
    has_priors = 0;

    num_robot_poses = problem_data.num_poses-has_priors; %remove the aux pose
    num_robots = problem_data.num_robots;
    poses_per_robot = num_robot_poses / num_robots;

    % skip the aux pose
    all_robot_start_idx = (problem_data.dim * has_priors)+1;
    per_robot_rot_block_span = 1:problem_data.dim * poses_per_robot;

    cur_robot_rot_block_start = (robot_idx-1) * problem_data.dim * poses_per_robot;
    all_robot_R_idxs = problem_data.all_R_idxs(all_robot_start_idx:end);
    R_idxs = all_robot_R_idxs(cur_robot_rot_block_start + per_robot_rot_block_span);
end