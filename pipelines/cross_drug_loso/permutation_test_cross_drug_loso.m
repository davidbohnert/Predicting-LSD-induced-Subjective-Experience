function output = permutation_test_cross_drug_loso( ...
        matsTrain, behavTrain, covarsTrain, trainGlobalRows, trainCohorts, ...
        heldoutTrainIdx, targetCohort, matsTestA, behavTestA, ...
        matsTestB, behavTestB, pThresholds, corrType, ...
        noIterations, useParallel, rngSettings)
%PERMUTATION_TEST_CROSS_DRUG_LOSO
% Subject-wise cross-drug LOSO CPM with one or more edge thresholds.
%
% The matching LSD subject is removed before edge selection and model
% fitting in every fold. If covariates are supplied, fold-specific partial
% correlation uses the sample-specific covariates selected by the caller:
% [study_2, study_3, sex, age, mean_FD] with study 1 as reference in the full
% sample, or [study_3, sex, age, mean_FD] with study 2 as reference for MEQ30.
% Adjustment is applied during edge selection. Test-set
% covariates are not used by the historical CPM prediction model.
%
% All requested thresholds share the same fixed-edge correlation cache.
% This changes computation order only: every threshold still receives its
% own edge masks, three CPM fits, predictions, and permutation inference.

if nargin < 15 || isempty(useParallel), useParallel = true; end
if nargin < 16 || isempty(rngSettings)
    rngSettings = struct('permutation_seed_offset',200, ...
        'loso_permutation_generator','twister');
end
assert(isfield(rngSettings, 'permutation_seed_offset') && ...
    isfield(rngSettings, 'loso_permutation_generator'), ...
    'CPM:RNG:MissingLOSOSettings', ...
    'LOSO permutation RNG settings are incomplete.');
repositoryDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repositoryDir, 'scripts', 'common'));

validateattributes(noIterations, {'numeric'}, ...
    {'scalar', 'integer', '>=', 2});
pThresholds = pThresholds(:)';
validateattributes(pThresholds, {'numeric'}, ...
    {'vector', 'real', 'finite', '>', 0, '<', 1});
assert(numel(unique(pThresholds)) == numel(pThresholds), ...
    'CPM:ExactLOSO:DuplicateThreshold', ...
    'Every edge threshold must be unique.');
corrType = char(validatestring(corrType, {'Pearson', 'Spearman'}));

behavTrain = behavTrain(:);
trainGlobalRows = trainGlobalRows(:);
trainCohorts = string(trainCohorts(:));
heldoutTrainIdx = heldoutTrainIdx(:);
behavTestA = behavTestA(:);
behavTestB = behavTestB(:);
targetCohort = string(targetCohort);

validate_inputs(matsTrain, behavTrain, covarsTrain, trainGlobalRows, ...
    trainCohorts, heldoutTrainIdx, targetCohort, ...
    matsTestA, behavTestA, matsTestB, behavTestB);

targetTrainIdx = find(trainCohorts == targetCohort);
assert(isequal(targetTrainIdx, heldoutTrainIdx), ...
    'CPM:ExactLOSO:TargetOrderMismatch', ...
    'Target-cohort LSD rows must match the two test-drug subject orders.');

nNodes = size(matsTrain, 1);
upperTriangle = triu(true(nNodes), 1);
edgeIndex = find(upperTriangle);
trainEdges = vectorize_unique_edges(matsTrain, edgeIndex);
testEdgesA = vectorize_unique_edges(matsTestA, edgeIndex);
testEdgesB = vectorize_unique_edges(matsTestB, edgeIndex);
nTest = numel(heldoutTrainIdx);
nEdges = numel(edgeIndex);
nThresholds = numel(pThresholds);

if isempty(covarsTrain), covarLabel = 'none';
else, covarLabel = sprintf('%d columns', size(covarsTrain, 2));
end
fprintf(['  Cached exact LOSO: cohort=%s, n_test=%d, n_train/fold=%d, ' ...
    'iterations=%d, thresholds=%s, covariates=%s\n'], ...
    targetCohort, nTest, size(trainEdges, 2) - 1, noIterations, ...
    mat2str(pThresholds), covarLabel);

% Column 1 is observed. Columns 2:end use the established subject-profile
% Configured subject-profile permutations; the publication convention is
% Twister with seeds 202:(noIterations+200).
permTrain = zeros(numel(behavTrain), noIterations);
permTestA = zeros(nTest, noIterations);
permTestB = zeros(nTest, noIterations);
permTrain(:, 1) = behavTrain;
permTestA(:, 1) = behavTestA;
permTestB(:, 1) = behavTestB;
for iteration = 2:noIterations
    [permTrain(:, iteration), permTestA(:, iteration), ...
        permTestB(:, iteration)] = permute_subject_profiles( ...
        behavTrain, behavTestA, behavTestB, trainCohorts, ...
        targetCohort, heldoutTrainIdx, ...
        iteration + rngSettings.permutation_seed_offset, ...
        rngSettings.loso_permutation_generator);
end

foldOutput = cell(nTest, 1);
if useParallel
    parfor fold = 1:nTest
        foldOutput{fold} = run_one_fold(fold, trainEdges, ...
            testEdgesA, testEdgesB, permTrain, covarsTrain, ...
            heldoutTrainIdx, pThresholds, corrType, noIterations);
    end
else
    for fold = 1:nTest
        foldOutput{fold} = run_one_fold(fold, trainEdges, ...
            testEdgesA, testEdgesB, permTrain, covarsTrain, ...
            heldoutTrainIdx, pThresholds, corrType, noIterations);
    end
end

predA = zeros(nTest, 3, nThresholds, noIterations);
predB = zeros(nTest, 3, nThresholds, noIterations);
baselinePred = zeros(nTest, noIterations);
posSelectionCount = zeros(nEdges, nThresholds, 'uint16');
negSelectionCount = zeros(nEdges, nThresholds, 'uint16');
posEdgeCountByFold = zeros(nTest, nThresholds);
negEdgeCountByFold = zeros(nTest, nThresholds);
for fold = 1:nTest
    predA(fold, :, :, :) = permute(foldOutput{fold}.pred_a, [4 2 3 1]);
    predB(fold, :, :, :) = permute(foldOutput{fold}.pred_b, [4 2 3 1]);
    baselinePred(fold, :) = foldOutput{fold}.baseline';
    posSelectionCount = posSelectionCount + ...
        uint16(foldOutput{fold}.observed_pos);
    negSelectionCount = negSelectionCount + ...
        uint16(foldOutput{fold}.observed_neg);
    posEdgeCountByFold(fold, :) = ...
        sum(foldOutput{fold}.observed_pos, 1);
    negEdgeCountByFold(fold, :) = ...
        sum(foldOutput{fold}.observed_neg, 1);
end

permRA = zeros(noIterations, 3, nThresholds);
permRB = zeros(noIterations, 3, nThresholds);
permRCombined = zeros(noIterations, 3, nThresholds);
permMSEA = zeros(noIterations, 3, nThresholds);
permMSEB = zeros(noIterations, 3, nThresholds);
permMSECombined = zeros(noIterations, 3, nThresholds);
for thresholdIdx = 1:nThresholds
    for iteration = 1:noIterations
        currentPredA = squeeze(predA(:, :, thresholdIdx, iteration));
        currentPredB = squeeze(predB(:, :, thresholdIdx, iteration));
        currentPredCombined = [currentPredA; currentPredB];
        currentTestCombined = [ ...
            permTestA(:, iteration); permTestB(:, iteration)];
        for modelIdx = 1:3
            permRA(iteration, modelIdx, thresholdIdx) = corr( ...
                currentPredA(:, modelIdx), permTestA(:, iteration), ...
                'Type', 'Pearson');
            permRB(iteration, modelIdx, thresholdIdx) = corr( ...
                currentPredB(:, modelIdx), permTestB(:, iteration), ...
                'Type', 'Pearson');
            permRCombined(iteration, modelIdx, thresholdIdx) = corr( ...
                currentPredCombined(:, modelIdx), currentTestCombined, ...
                'Type', 'Pearson');
            permMSEA(iteration, modelIdx, thresholdIdx) = mean(( ...
                permTestA(:, iteration) - currentPredA(:, modelIdx)).^2);
            permMSEB(iteration, modelIdx, thresholdIdx) = mean(( ...
                permTestB(:, iteration) - currentPredB(:, modelIdx)).^2);
            permMSECombined(iteration, modelIdx, thresholdIdx) = mean(( ...
                currentTestCombined - ...
                currentPredCombined(:, modelIdx)).^2);
        end
    end
end

assert(all(isfinite(permRA(:))) && all(isfinite(permRB(:))) && ...
    all(isfinite(permRCombined(:))) && ...
    all(isfinite(permMSEA(:))) && all(isfinite(permMSEB(:))) && ...
    all(isfinite(permMSECombined(:))), ...
    'CPM:ExactLOSO:NonfinitePerformance', ...
    'The LOSO pipeline produced a nonfinite performance statistic.');

output.thresholds = pThresholds;
observedA = cell(nThresholds, 1);
observedB = cell(nThresholds, 1);
observedCombined = cell(nThresholds, 1);
for thresholdIdx = 1:nThresholds
    currentPredA = squeeze(predA(:, :, thresholdIdx, 1));
    currentPredB = squeeze(predB(:, :, thresholdIdx, 1));
    observedA{thresholdIdx} = make_observed_stats( ...
        currentPredA, baselinePred(:, 1), behavTestA, ...
        permRA(1, :, thresholdIdx), permMSEA(1, :, thresholdIdx));
    observedB{thresholdIdx} = make_observed_stats( ...
        currentPredB, baselinePred(:, 1), behavTestB, ...
        permRB(1, :, thresholdIdx), permMSEB(1, :, thresholdIdx));
    observedCombined{thresholdIdx} = make_observed_stats( ...
        [currentPredA; currentPredB], ...
        [baselinePred(:, 1); baselinePred(:, 1)], ...
        [behavTestA; behavTestB], ...
        permRCombined(1, :, thresholdIdx), ...
        permMSECombined(1, :, thresholdIdx));
    currentDrugA = add_permutation_pvalues( ...
        observedA{thresholdIdx}, permRA(:, :, thresholdIdx), ...
        permMSEA(:, :, thresholdIdx), noIterations);
    currentDrugB = add_permutation_pvalues( ...
        observedB{thresholdIdx}, permRB(:, :, thresholdIdx), ...
        permMSEB(:, :, thresholdIdx), noIterations);
    currentCombined = add_permutation_pvalues( ...
        observedCombined{thresholdIdx}, ...
        permRCombined(:, :, thresholdIdx), ...
        permMSECombined(:, :, thresholdIdx), noIterations);
    if thresholdIdx == 1
        output.drug_a = currentDrugA;
        output.drug_b = currentDrugB;
        output.combined = currentCombined;
    else
        output.drug_a(thresholdIdx) = currentDrugA;
        output.drug_b(thresholdIdx) = currentDrugB;
        output.combined(thresholdIdx) = currentCombined;
    end
end

audit = struct;
audit.train_global_rows = trainGlobalRows;
audit.heldout_train_indices = heldoutTrainIdx;
audit.heldout_global_rows = trainGlobalRows(heldoutTrainIdx);
audit.train_cohorts = trainCohorts;
audit.target_cohort = targetCohort;
audit.no_nodes = nNodes;
audit.edge_index = edgeIndex;
audit.p_thresholds = pThresholds;
audit.correlation_type = corrType;
audit.performance_correlation = 'Pearson';
audit.no_iterations = noIterations;
audit.permutation_seeds = ((2:noIterations)' + ...
    rngSettings.permutation_seed_offset);
audit.permutation_generator = rngSettings.loso_permutation_generator;
audit.n_train_per_fold = size(trainEdges, 2) - 1;
audit.n_test_single_drug = nTest;
audit.n_test_combined = 2 * nTest;
audit.pos_selection_count = posSelectionCount;
audit.neg_selection_count = negSelectionCount;
audit.pos_selection_frequency = double(posSelectionCount) / nTest;
audit.neg_selection_frequency = double(negSelectionCount) / nTest;
audit.pos_edge_count_by_fold = posEdgeCountByFold;
audit.neg_edge_count_by_fold = negEdgeCountByFold;
audit.observed_a = observedA;
audit.observed_b = observedB;
audit.observed_combined = observedCombined;
audit.null_r_a = permRA(2:end, :, :);
audit.null_r_b = permRB(2:end, :, :);
audit.null_r_combined = permRCombined(2:end, :, :);
audit.null_mse_a = permMSEA(2:end, :, :);
audit.null_mse_b = permMSEB(2:end, :, :);
audit.null_mse_combined = permMSECombined(2:end, :, :);
audit.covariates_used = ~isempty(covarsTrain);
audit.covariate_column_count = size(covarsTrain, 2);
output.audit = audit;
end


function result = run_one_fold(fold, trainEdges, testEdgesA, ...
        testEdgesB, permTrain, covarsTrain, heldoutTrainIdx, ...
        pThresholds, corrType, noIterations)
keepTrain = true(size(trainEdges, 2), 1);
keepTrain(heldoutTrainIdx(fold)) = false;
foldEdges = trainEdges(:, keepTrain);
if isempty(covarsTrain), foldCovars = [];
else, foldCovars = covarsTrain(keepTrain, :);
end
cache = cpm_prepare_correlation_cache(foldEdges', foldCovars, corrType);
nThresholds = numel(pThresholds);
nEdges = size(trainEdges, 1);
predA = zeros(noIterations, 3, nThresholds);
predB = zeros(noIterations, 3, nThresholds);
baseline = zeros(noIterations, 1);
observedPos = false(nEdges, nThresholds);
observedNeg = false(nEdges, nThresholds);

for iteration = 1:noIterations
    foldBehav = permTrain(keepTrain, iteration);
    [~, posMasks, negMasks] = cpm_cached_edge_selection( ...
        cache, foldBehav, pThresholds);
    for thresholdIdx = 1:nThresholds
        posMask = posMasks{thresholdIdx};
        negMask = negMasks{thresholdIdx};
        trainSumPos = sum(foldEdges(posMask, :), 1)';
        trainSumNeg = sum(foldEdges(negMask, :), 1)';
        [bComb, bPos, bNeg] = cpm_fit_three_models_min_norm( ...
            foldBehav, trainSumPos, trainSumNeg);
        testASumPos = sum(testEdgesA(posMask, fold), 1);
        testASumNeg = sum(testEdgesA(negMask, fold), 1);
        testBSumPos = sum(testEdgesB(posMask, fold), 1);
        testBSumNeg = sum(testEdgesB(negMask, fold), 1);
        predA(iteration, :, thresholdIdx) = [ ...
            bComb(1) * testASumPos + bComb(2) * testASumNeg + bComb(3), ...
            bPos(1) * testASumPos + bPos(2), ...
            bNeg(1) * testASumNeg + bNeg(2)];
        predB(iteration, :, thresholdIdx) = [ ...
            bComb(1) * testBSumPos + bComb(2) * testBSumNeg + bComb(3), ...
            bPos(1) * testBSumPos + bPos(2), ...
            bNeg(1) * testBSumNeg + bNeg(2)];
        if iteration == 1
            observedPos(:, thresholdIdx) = posMask;
            observedNeg(:, thresholdIdx) = negMask;
        end
    end
    baseline(iteration) = mean(foldBehav);
end
result = struct('pred_a', predA, 'pred_b', predB, ...
    'baseline', baseline, 'observed_pos', observedPos, ...
    'observed_neg', observedNeg);
if mod(fold, 5) == 0 || fold == numel(heldoutTrainIdx)
    fprintf('  ...LOSO fold %d / %d complete\n', ...
        fold, numel(heldoutTrainIdx));
end
end


function observed = make_observed_stats(predictions, baselinePred, ...
        behavior, observedR, observedMSE)
observed = struct;
observed.r = reshape(observedR, 1, []);
observed.mse = reshape(observedMSE, 1, []);
observed.baseline_mse = mean((behavior - baselinePred).^2);
observed.pred_comb = predictions(:, 1);
observed.pred_pos = predictions(:, 2);
observed.pred_neg = predictions(:, 3);
observed.baseline_pred = baselinePred;
observed.observed = behavior;
end


function result = add_permutation_pvalues(observed, permR, permMSE, nIter)
result = observed;
result.p_r = [sum(permR(:, 1) >= observed.r(1)), ...
              sum(permR(:, 2) >= observed.r(2)), ...
              sum(permR(:, 3) >= observed.r(3))] / nIter;
result.p_mse = [sum(permMSE(:, 1) <= observed.mse(1)), ...
                sum(permMSE(:, 2) <= observed.mse(2)), ...
                sum(permMSE(:, 3) <= observed.mse(3))] / nIter;
end


function [permutedTrain, permutedTestA, permutedTestB] = ...
        permute_subject_profiles(trainBehav, testBehavA, testBehavB, ...
        trainCohorts, targetCohort, heldoutTrainIdx, seed, generator)
rng(seed, generator);
permutedTrain = trainBehav;
permutedTestA = testBehavA;
permutedTestB = testBehavB;
cohortOrder = ["5HT", "LAM", "LPM"];
for cohort = cohortOrder
    cohortIdx = find(trainCohorts == cohort);
    if isempty(cohortIdx), continue; end
    sourceOrder = randperm(numel(cohortIdx));
    permutedTrain(cohortIdx) = trainBehav(cohortIdx(sourceOrder));
    if cohort == targetCohort
        assert(isequal(cohortIdx, heldoutTrainIdx), ...
            'CPM:ExactLOSO:PermutationTargetOrderMismatch', ...
            'Target behavior profiles are not aligned to LSD rows.');
        permutedTestA = testBehavA(sourceOrder);
        permutedTestB = testBehavB(sourceOrder);
    end
end
end


function edges = vectorize_unique_edges(mats, edgeIndex)
flattened = reshape(mats, [], size(mats, 3));
edges = flattened(edgeIndex, :);
end


function validate_inputs(matsTrain, behavTrain, covarsTrain, ...
        trainGlobalRows, trainCohorts, heldoutTrainIdx, targetCohort, ...
        matsTestA, behavTestA, matsTestB, behavTestB)
values = {matsTrain, matsTestA, matsTestB};
for idx = 1:numel(values)
    mats = values{idx};
    assert(isnumeric(mats) && ndims(mats) == 3 && ...
        size(mats, 1) == size(mats, 2) && all(isfinite(mats(:))), ...
        'CPM:ExactLOSO:InvalidConnectome', ...
        'Every connectome input must be a finite square 3-D array.');
end
assert(size(matsTrain, 1) == size(matsTestA, 1) && ...
    size(matsTrain, 1) == size(matsTestB, 1), ...
    'CPM:ExactLOSO:NodeCountMismatch', ...
    'Training and testing connectomes have different node counts.');
assert(size(matsTrain, 3) == numel(behavTrain) && ...
    numel(trainGlobalRows) == numel(behavTrain) && ...
    numel(trainCohorts) == numel(behavTrain), ...
    'CPM:ExactLOSO:TrainingAlignment', ...
    'Training connectomes, behavior, and metadata are not aligned.');
assert(numel(heldoutTrainIdx) == numel(behavTestA) && ...
    numel(heldoutTrainIdx) == numel(behavTestB) && ...
    size(matsTestA, 3) == numel(behavTestA) && ...
    size(matsTestB, 3) == numel(behavTestB), ...
    'CPM:ExactLOSO:TestingAlignment', ...
    'Held-out mappings and test inputs are not aligned.');
assert(all(isfinite(behavTrain)) && all(isfinite(behavTestA)) && ...
    all(isfinite(behavTestB)), ...
    'CPM:ExactLOSO:NonfiniteBehavior', ...
    'Behavior inputs contain NaN or Inf.');
if ~isempty(covarsTrain)
    assert(isnumeric(covarsTrain) && ...
        size(covarsTrain, 1) == numel(behavTrain) && ...
        all(isfinite(covarsTrain(:))), ...
        'CPM:ExactLOSO:InvalidCovariates', ...
        'Covariates must be finite and aligned to LSD training rows.');
end
assert(all(trainCohorts(heldoutTrainIdx) == targetCohort), ...
    'CPM:ExactLOSO:HeldoutCohort', ...
    'Held-out LSD rows do not belong to the target cohort.');
end
