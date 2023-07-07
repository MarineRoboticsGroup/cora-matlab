% get the initialization
function lifted_init = lift_init_point(problem, X, lifted_manifold, add_noise, saddle_escape_dir, ~, tolgradnorm)

    % require that all arguments are given
    if nargin < 6
        error('lift_init_point requires at least %d arguments but was given %d ', expected_nargin, nargin);
    end

    % if add_noise then saddle_escape_dir must be []
    if add_noise && ~isempty(saddle_escape_dir)
        error('lift_init_point: if add_noise is true, then saddle_escape_dir must be []');
    end

    base_dim = lifted_manifold.base_dim;
    lift_dim = lifted_manifold.lifted_dim;

    % make sure base_dim <= size(X, 2) <= lift_dim
    assert(base_dim <= size(X, 2) && size(X, 2) <= lift_dim, ...
        'lift_init_point: X must have size(X, 2) in [base_dim, lifted_dim] but has shape %s', mat2str(size(X)));


    % pad X with zeros to make it the right shape
    lifted_init = lifted_manifold.zeros();
    lifted_init(:, 1:size(X, 2)) = X;

    % if saddle_escape_dir is given, then add_noise must be false and we are only
    % incrementing the size of X by one
    if ~isempty(saddle_escape_dir)
        assert(size(saddle_escape_dir, 2) == 1, 'lift_init_point: if saddle_escape_dir is given, then size(saddle_escape_dir, 2) must be 1');
        assert(size(X, 2) + 1 == lift_dim, 'lift_init_point: if saddle_escape_dir is given, then size(X, 2) + 1 must be lift_dim');
    elseif add_noise
        % if desired, add some noise to the lifted point. This is useful becausei
        % if rank(lifted_init) < lifted_dim, then it is likely (due to the
        % fundamentals of manifold optimization) that we will not be able to
        % estimate a higher rank point
        saddle_escape_dir = rand(size(lifted_init));
    end

    if add_noise || ~isempty(saddle_escape_dir)
        % if we have a saddle_escape_dir, treat it as a tangent vector and retract
        % lifted_init = lifted_manifold.retr(lifted_init, saddle_escape_dir, 1e-2);
        [lifted_init, search_success] = retract_along_direction_with_linesearch(problem, tolgradnorm, lifted_init, saddle_escape_dir);
    end

    % randomly apply a rotation to the lifted point to fill the zero entries
    % but keeping the point the same (up to SO gauge symmetry)
    rot = randrot(lift_dim);
    lifted_init = lifted_init * rot;
end

