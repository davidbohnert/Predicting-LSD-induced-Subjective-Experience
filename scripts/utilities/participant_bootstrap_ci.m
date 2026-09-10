function interval = participant_bootstrap_ci(observed, predicted, ...
        noBootstraps, seed, profileIds)
%PARTICIPANT_BOOTSTRAP_CI Participant-level percentile and BCa intervals.
%
% PROFILEIDS groups repeated sessions from the same participant. Omitting it
% treats each row as an independent participant. All rows belonging to a
% sampled profile are retained together.

if nargin < 3 || isempty(noBootstraps), noBootstraps = 10000; end
if nargin < 5 || isempty(profileIds), profileIds = (1:numel(observed))'; end
observed = observed(:); predicted = predicted(:); profileIds = profileIds(:);
assert(numel(observed) == numel(predicted) && ...
    numel(observed) == numel(profileIds), 'CPM:Bootstrap:InputSize', ...
    'Observed, predicted, and profile identifiers must align.');

profiles = unique(profileIds, 'stable');
nProfiles = numel(profiles);
if nargin < 4 || isempty(seed), seed = 1300000 + nProfiles; end
stream = RandStream('mt19937ar', 'Seed', seed);
% Preserve the resampling convention: MATLAB fills this matrix
% column-major, and each row is one participant bootstrap replicate.
draws = randi(stream, nProfiles, noBootstraps, nProfiles);
replicates = nan(noBootstraps, 1);
for iteration = 1:noBootstraps
    sampled = profiles(draws(iteration, :));
    blocks = arrayfun(@(id) find(profileIds == id), sampled, ...
        'UniformOutput', false);
    indices = vertcat(blocks{:});
    replicates(iteration) = corr(observed(indices), predicted(indices), ...
        'Type', 'Pearson');
end

point = corr(observed, predicted, 'Type', 'Pearson');
percentile = prctile(replicates, [2.5 97.5]);
jackknife = nan(nProfiles, 1);
for index = 1:nProfiles
    keep = profileIds ~= profiles(index);
    jackknife(index) = corr(observed(keep), predicted(keep), 'Type', 'Pearson');
end
z0 = norminv(mean(replicates < point));
centered = mean(jackknife) - jackknife;
acceleration = sum(centered.^3) / ...
    (6 * max(sum(centered.^2)^(3/2), eps));
z = norminv([0.025 0.975]);
adjusted = normcdf(z0 + (z0 + z) ./ (1 - acceleration * (z0 + z)));
bca = quantile(replicates, adjusted);
interval = struct('point',point, 'bca',bca, 'percentile',percentile, ...
    'replicates',replicates, 'profiles',nProfiles, ...
    'no_bootstraps',noBootstraps, 'seed',seed, 'draws',draws);
end
