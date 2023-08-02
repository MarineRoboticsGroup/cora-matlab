function rounded_soln = round_solution(X, problem, verbosity, dim)
    base_dim = problem.dim;
    if nargin < 4
        dim = problem.dim;
    end
    num_poses = problem.num_poses;

    % get U, S, Vh from the SVD (rounding is thin SVD)
    [U, S, V] = svd(X, "econ");
    Vh = V';

    % this check makes sure that we're multiplying the right components
    lifted_dim = size(X, 2);
    expected_Vh_size = [lifted_dim, lifted_dim];
    assert(isequal(size(Vh), expected_Vh_size), ...
        sprintf('The shape of Vh is %d x %d but should be %d x %d, the solution should be transposed\n', ...
        size(Vh), expected_Vh_size) ...
    );

    % round the solution via thin SVD
    rounded_soln = (U(:, 1:dim) * S(1:dim, 1:dim))';

    % simultaneously round the rotations in the solution (via SVD) and check if
    % they are left-handed, in which case they need to be flipped
    determinants = zeros(1, num_poses);
    for i = 0:num_poses-1
        rot_i_idxs = problem.all_R_idxs((i*base_dim)+1:((i+1)*base_dim)-1);
        rot_i = rounded_soln(:, rot_i_idxs);
        [Urot, S, Vrot] = svd(rot_i, "econ");
        Vhrot = Vrot';
        rot_i = Urot * Vhrot;
        rounded_soln(:, rot_i_idxs) = rot_i;
        determinants(i+1) = det(S);
    end

    num_good_rotations = sum(abs(determinants - 1) < 1e-2);
    if verbosity > 0
        fprintf('Number of good rotations: %d / %d\n',...
            num_good_rotations, num_poses);
    end

    if num_good_rotations == 0
        if verbosity > 0
            fprintf('All rotations were left handed, flipping them!\n');
        end
        reflector = diag([ones(1, dim - 1), -1]);
        rounded_soln = reflector * rounded_soln;
    elseif num_good_rotations < num_poses / 2
        if verbosity > 0
            fprintf('Majority of rotations were left handed, flipping them!\n');
        end
        reflector = diag([ones(1, dim - 1), -1]);
        rounded_soln = reflector * rounded_soln;
    end

    %% also round the distance variables to the nearest unit vector
    dist_idxs = problem.all_d_idxs;
    for i = 1:length(dist_idxs)
        dist_i = rounded_soln(:, dist_idxs(i));
        rounded_soln(:, dist_idxs(i)) = dist_i / norm(dist_i);
    end

end

