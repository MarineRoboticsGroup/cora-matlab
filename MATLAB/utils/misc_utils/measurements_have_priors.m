function has_priors = measurements_have_priors(measurements)
    % check if measurements.pose_priors or measurements.landmark_priors exist and are not empty
    has_pose_priors = isfield(measurements, 'pose_priors') && ~isempty(measurements.pose_priors.mean);
    has_landmark_priors = isfield(measurements, 'landmark_priors') && ~isempty(measurements.landmark_priors.mean);
    has_priors = has_pose_priors || has_landmark_priors;
end