function t_idxs = get_robot_t_idxs(problem_data, t_idxs, robot_idx)

    has_priors = (problem_data.num_pose_priors + problem_data.num_beacon_priors) > 0;

    num_robot_poses = problem_data.num_poses-has_priors; %remove the aux pose
    num_robots = problem_data.num_robots;
    poses_per_robot = num_robot_poses / num_robots;

    % skip the aux pose
    all_robot_t_idxs = t_idxs(has_priors+1:end);
    t_idxs = all_robot_t_idxs(poses_per_robot*(robot_idx-1)+1:poses_per_robot*robot_idx);

end