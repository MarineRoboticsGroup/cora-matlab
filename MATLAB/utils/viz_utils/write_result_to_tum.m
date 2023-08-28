function write_result_to_tum(X, problem_data, save_dir)

    if (problem_data.use_marginalized)
        Xfull = extract_translations_from_marginalized_solution(X, problem_data);
    else
        Xfull = X;
    end

    % if X is tall and skinny, transpose it
    if size(Xfull, 1) > size(Xfull, 2)
        Xfull = Xfull';
    end

    warning("Function does not currently correctly handle priors!!")
    % has_priors = (problem_data.num_pose_priors + problem_data.num_landmark_priors) > 0;
    has_priors = 0;

    dim = problem_data.dim;

    num_robot_poses = problem_data.num_poses-has_priors; %remove the aux pose
    num_robots = problem_data.num_robots;
    assert(mod(num_robot_poses, num_robots) == 0, 'Number of robot poses (%f) must be divisible by number of robots (%f)', num_robot_poses, num_robots);
    poses_per_robot = num_robot_poses / num_robots;
    num_rot_idxs_per_robot = poses_per_robot * dim;

    if num_robots > 0

        % make a list of the paths to save to, numbered by robot
        save_fpaths = save_dir + "/cora_" + string(1:num_robots) + ".tum";
        for robot_idx = 1:num_robots
            % open the file to write to
            fid = fopen(save_fpaths(robot_idx), 'w');
            if fid == -1
                error('Cannot open file: %s', save_fpaths(robot_idx));
            end

            cur_robot_trans_idxs = get_robot_t_idxs(problem_data, robot_idx);
            cur_robot_rot_idxs = get_robot_R_idxs(problem_data, robot_idx);

            trans_i = Xfull(:, cur_robot_trans_idxs);
            rotations_i = Xfull(:, cur_robot_rot_idxs);

            % if problem_data.timestamps exists, use it
            % otherwise, use the default timestamps
            if isfield(problem_data, 'timestamps')
                timestamp_start = (robot_idx - 1)*poses_per_robot+1;
                timestamp_end = timestamp_start + poses_per_robot - 1;
                robot_i_timestamps = problem_data.timestamps(timestamp_start:timestamp_end);
            else
                robot_i_timestamps = 0:poses_per_robot-1;
            end
            assert(length(robot_i_timestamps) == poses_per_robot);

            % write the poses to the file in the format:
            % timestamp x y z qx qy qz qw
            for fpath_line_idx = 1:length(trans_i)
                timestamp = robot_i_timestamps(fpath_line_idx);
                trans_i_j = trans_i(:, fpath_line_idx);
                x = trans_i_j(1);
                y = trans_i_j(2);
                if(dim == 2)
                    z = 0;
                elseif dim == 3
                    z = trans_i_j(3);
                end

                rot_j_start_idx = (fpath_line_idx-1)*dim+1;
                rot_j_idxs = rot_j_start_idx:rot_j_start_idx+dim-1;
                rot_i_j = rotations_i(:, rot_j_idxs);

                % convert the rotation matrix to a quaternion
                % if rotation is 2x2, pad it with a row and column of zeros
                if dim == 2
                    rot_i_j_quat = rotm2quat([rot_i_j, zeros(2,1); zeros(1,3)]);
                elseif dim == 3
                    rot_i_j_quat = rotm2quat(rot_i_j);
                else
                    error('Invalid dimension: %d', dim);
                end
                qw = rot_i_j_quat(1);
                qx = rot_i_j_quat(2);
                qy = rot_i_j_quat(3);
                qz = rot_i_j_quat(4);

                % write the line to the file
                fprintf(fid, '%f %f %f %f %f %f %f %f\n', timestamp, x, y, z, qx, qy, qz, qw);

            end

            % close the file
            fclose(fid);
            % fprintf('Saved robot poses to %s\n', save_fpaths(robot_idx));
        end
    end



end
