function test_optimized_equivalence
%TEST_OPTIMIZED_EQUIVALENCE Compare readable and cached internal CPM cores.
testDir = fileparts(mfilename('fullpath'));
repositoryRoot = fileparts(testDir);
addpath(fullfile(repositoryRoot, 'tests', 'reference'), ...
    fullfile(repositoryRoot, 'scripts', 'common'), ...
    fullfile(repositoryRoot, 'pipelines', 'internal_validation'));

rng(901, 'twister');
nNodes = 10; nSubjects = 24;
mats = zeros(nNodes, nNodes, nSubjects);
for subject = 1:nSubjects
    raw = randn(nNodes);
    current = (raw + raw') / 2;
    current(1:nNodes+1:end) = 0;
    mats(:,:,subject) = current;
end
behavior = randn(nSubjects, 1);
adjustmentSets = {[], [randn(nSubjects,1), double((1:nSubjects)'>12)]};
thresholds = [0.10 0.20];
for corrType = {'Pearson','Spearman'}
    for adjustmentIndex = 1:numel(adjustmentSets)
        covars = adjustmentSets{adjustmentIndex};
        textPath = [tempname '.txt'];
        cleanup = onCleanup(@() delete_if_present(textPath));
        [rows, details] = permutation_test_cv_optimized('synthetic', ...
            'synthetic_behavior', mats, behavior, 6, [0.10 0.20], ...
            covars, textPath, {}, 2, corrType{1}, 0);
        for thresholdIndex = 1:2
            rng(123, 'twister');
            partition = cvpartition(nSubjects, 'KFold', 6);
            [rCombined,rPositive,rNegative,mseCombined,msePositive,mseNegative] = ...
                CPM_core_function_reference(mats, behavior, partition, ...
                    thresholds(thresholdIndex), covars, corrType{1});
            optimized = cell2mat(rows(thresholdIndex,[4 6 8 10 12 14]));
            reference = [rCombined,rPositive,rNegative, ...
                mseCombined,msePositive,mseNegative];
            assert(max(abs(optimized-reference)) < 1e-10, ...
                'CPM:Test:InternalEquivalence', ...
                'Optimized and readable internal CPM results differ.');
        end
        assert(details.uses_unique_edges && ...
            isequal(size(details.observed_predictions.positive), [nSubjects 2]));
        clear cleanup
    end
end

customRng = struct('cv_seed',321, 'cv_generator','twister', ...
    'permutation_seed_offset',700, ...
    'internal_permutation_generator','threefry');
textPath = [tempname '.txt'];
cleanup = onCleanup(@() delete_if_present(textPath));
[~, customDetails] = permutation_test_cv_optimized('synthetic', ...
    'synthetic_behavior', mats, behavior, 6, 0.10, [], textPath, {}, ...
    2, 'Pearson', 0, customRng);
rng(customRng.cv_seed, customRng.cv_generator);
expectedPartition = cvpartition(nSubjects, 'KFold', 6);
for fold = 1:6
    assert(isequal(test(customDetails.cv_partition, fold), ...
        test(expectedPartition, fold)), ...
        'Configured CV seed did not reach the optimized engine.');
end
expectedPermutations = cpm_behavior_permutations(behavior, 2, customRng);
assert(isequal(customDetails.perm_behav, expectedPermutations), ...
    'Configured permutation settings did not reach the optimized engine.');
clear cleanup
fprintf('Pearson/Spearman, adjusted/unadjusted, multi-threshold equivalence passed.\n');
end

function delete_if_present(pathName)
if isfile(pathName), delete(pathName); end
end
