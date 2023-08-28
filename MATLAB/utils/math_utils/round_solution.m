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
        offset = i * base_dim;
        rot_i_idxs = problem.all_R_idxs(offset + (1:base_dim));
        rot_i = rounded_soln(:, rot_i_idxs);
        [Urot, S, Vrot] = svd(rot_i, "econ");
        Vhrot = Vrot';
        rot_i = Urot * Vhrot;

        % rot_i should be dim x dim
        assert(isequal(size(rot_i), [dim, dim]), ...
            sprintf('The shape of rot_i is %d x %d but should be %d x %d\n', ...
            size(rot_i), [dim, dim]) ...
        );

        rounded_soln(:, rot_i_idxs) = rot_i;
        determinants(i+1) = det(rot_i);
        % determinant should be 1 or -1, but we allow for some numerical error
        assert((abs(determinants(i+1)) - 1) < 1e-6, 'Determinant of rotation %d is %f, should be 1 or -1\n', ...
            i, determinants(i+1));
    end

    num_good_rotations = sum(abs(determinants - 1) < 1e-6);
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

        % make sure dist_i is a vector of length dim
        assert(isequal(size(dist_i), [dim, 1]), ...
            sprintf('Distance %d is not a vector of length %d\n', i, dim));

        % check that the distance is a unit vector
        assert(abs(norm(rounded_soln(:, dist_idxs(i))) - 1) < 1e-6, ...
            sprintf('Distance %d is not a unit vector\n', i));
    end

end

