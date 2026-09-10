function [result, details] = permutation_test_cv_study_loso_positive( ...
        matrixName, behaviorName, mats, behavior, studyLabels, covars, ...
        threshold, noIterations, corrType, maxWorkers, rngSettings)
%PERMUTATION_TEST_CV_STUDY_LOSO_POSITIVE Study-held-out, adjusted positive CPM.
behavior = behavior(:);
studyLabels = studyLabels(:);
n = size(mats, 3);
if nargin < 10 || isempty(maxWorkers)
    maxWorkers = 0;
end
if nargin < 11 || isempty(rngSettings)
    rngSettings = struct('permutation_seed_offset',200, ...
        'internal_permutation_generator','threefry');
end
assert(numel(behavior)==n && numel(studyLabels)==n && size(covars,1)==n, ...
    'CPM:PrimaryStudyLOSO:Alignment', ...
    'Connectomes, behavior, studies, and covariates must align.');
assert(all(isfinite(mats),'all') && all(isfinite(behavior)) && all(isfinite(covars),'all'), ...
    'CPM:PrimaryStudyLOSO:Finite', ...
    'All analysis inputs must be finite.');
studyIds = unique(studyLabels, 'stable');
assert(numel(studyIds) >= 2, 'CPM:PrimaryStudyLOSO:Folds', ...
    'At least two studies are required.');
fprintf(['Running study LOSO: %s vs %s; %s selection; ' ...
    'p<%.3f; N=%d\n'], ...
    matrixName, behaviorName, corrType, threshold, n);
permutedBehavior = cpm_behavior_permutations( ...
    behavior, noIterations, rngSettings);
layout = cpm_prepare_edge_layout(mats, corrType);
foldOutput = cell(numel(studyIds),1);
parfor (fold = 1:numel(studyIds), maxWorkers)
    trainIndex = find(studyLabels ~= studyIds(fold));
    testIndex = find(studyLabels == studyIds(fold));
    [foldCovars, foldNames, foldReference] = cpm_select_covariates( ...
        covars(trainIndex,:), studyLabels(trainIndex));
    cache = cpm_prepare_correlation_cache( ...
        layout.correlation_edges(:,trainIndex)', foldCovars, corrType);
    foldOutput{fold} = run_fold( ...
        layout.strength_edges(:,trainIndex), ...
        layout.strength_edges(:,testIndex), layout.strength_divisor, ...
        trainIndex, testIndex, permutedBehavior, threshold, cache);
    foldOutput{fold}.covariate_columns = foldNames;
    foldOutput{fold}.study_reference = foldReference;
    fprintf('  held-out study %g: %d test subjects\n', studyIds(fold), numel(testIndex));
end
predicted = zeros(n, noIterations);
for fold = 1:numel(studyIds)
    predicted(foldOutput{fold}.testIndex,:) = foldOutput{fold}.predicted;
end
permutationR = zeros(noIterations,1);
permutationMse = zeros(noIterations,1);
for iteration = 1:noIterations
    current = permutedBehavior(:,iteration);
    permutationR(iteration) = corr(predicted(:,iteration), current);
    permutationMse(iteration) = mean((current-predicted(:,iteration)).^2);
end
assert(all(isfinite(permutationR)) && all(isfinite(permutationMse)), ...
    'CPM:PrimaryStudyLOSO:Performance', ...
    'Nonfinite performance statistic.');
result = struct('r',permutationR(1), ...
    'p',sum(permutationR>=permutationR(1))/noIterations, ...
    'mse',permutationMse(1));
details = struct('permutation_r',permutationR,'permutation_mse',permutationMse, ...
    'study_ids',studyIds,'study_labels',studyLabels,'corr_type',corrType, ...
    'threshold',threshold,'uses_unique_edges',layout.uses_unique_edges, ...
    'fold_details',{foldOutput});
end

function output = run_fold(trainEdges, testEdges, divisor, trainIndex, testIndex, permutedBehavior, threshold, cache)
iterations = size(permutedBehavior,2);
predicted = zeros(numel(testIndex),iterations);
for iteration = 1:iterations
    foldBehavior = permutedBehavior(trainIndex,iteration);
    [~, posMasks] = cpm_cached_edge_selection( ...
        cache, foldBehavior, threshold);
    posMask = posMasks{1};
    trainStrength = sum(trainEdges(posMask,:),1)' / divisor;
    testStrength = sum(testEdges(posMask,:),1)' / divisor;
    [~, bPos] = cpm_fit_three_models( ...
        foldBehavior, trainStrength, zeros(size(trainStrength)));
    predicted(:,iteration) = bPos(1) * testStrength + bPos(2);
end
output = struct('testIndex',testIndex,'predicted',predicted);
end
