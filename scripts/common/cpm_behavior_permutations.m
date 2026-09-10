function permuted_behavior = cpm_behavior_permutations( ...
        behavior, no_iterations, rng_settings)
% CPM_BEHAVIOR_PERMUTATIONS
% Reproduce the observed-plus-shuffled vectors used by the original tests.
%
% The submitted wrappers generated each shuffle inside their PARFOR loops.
% Caching is easier when the same vectors are prepared once on the client,
% so this helper changes where—not how—the permutations are constructed.
% The original process-worker calls used Threefry with seed iteration+200.
% RNG_SETTINGS makes those fixed conventions explicit and auditable.

if nargin < 3 || isempty(rng_settings)
    rng_settings = struct('permutation_seed_offset', 200, ...
        'internal_permutation_generator', 'threefry');
end
assert(isfield(rng_settings, 'permutation_seed_offset') && ...
    isfield(rng_settings, 'internal_permutation_generator'), ...
    'CPM:RNG:MissingInternalSettings', ...
    'Internal permutation RNG settings are incomplete.');

behavior = behavior(:);
no_subjects = numel(behavior);
permuted_behavior = zeros( ...
    no_subjects, no_iterations, 'like', behavior);
permuted_behavior(:, 1) = behavior;

for iteration = 2:no_iterations
    stream = RandStream(rng_settings.internal_permutation_generator, ...
        'Seed', iteration + rng_settings.permutation_seed_offset);
    order = randperm(stream, no_subjects);
    permuted_behavior(:, iteration) = behavior(order);
end
end
