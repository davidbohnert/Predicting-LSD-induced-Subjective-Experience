function plan = learning_curve_plan(planPath, standardLabels, meqLabels, ...
        standardSampleSizes, meqSampleSizes, noRepetitions, rngSettings)
%LEARNING_CURVE_PLAN Reproducible nested study-stratified samples.
%
% Within a repetition, every smaller sample is contained in the next larger
% sample. The same samples and CV seeds are reused across connectome forms,
% behavior forms, and covariate settings so condition comparisons are paired.

standardLabels = standardLabels(:);
meqLabels = meqLabels(:);
if nargin < 4 || isempty(standardSampleSizes)
    standardSampleSizes = [20 30 40 50 60 67];
end
if nargin < 5 || isempty(meqSampleSizes)
    meqSampleSizes = [20 30 40 48];
end
if nargin < 6 || isempty(noRepetitions), noRepetitions = 100; end
if nargin < 7 || isempty(rngSettings)
    rngSettings = struct('cv_seed',123, 'cv_generator','twister');
end
assert(numel(standardLabels) == 67 && ...
    isequal(unique(standardLabels)', [1 2 3]), ...
    'CPM:LearningCurve:StandardLabels', ...
    'Expected 67 standard-cohort labels spanning studies 1, 2, and 3.');
assert(numel(meqLabels) == 48 && ...
    isequal(unique(meqLabels)', [2 3]), ...
    'CPM:LearningCurve:MEQLabels', ...
    'Expected 48 MEQ30-cohort labels spanning studies 2 and 3.');

if isfile(planPath)
    loaded = load(planPath, 'plan');
    plan = loaded.plan;
    verify_plan(plan, standardLabels, meqLabels, standardSampleSizes, ...
        meqSampleSizes, noRepetitions, rngSettings);
    return
end

plan.version = 3;
plan.no_repetitions = noRepetitions;
plan.cv_seeds = rngSettings.cv_seed + (0:(noRepetitions - 1));
plan.cv_generator = rngSettings.cv_generator;
plan.sampling_generator = 'twister';
plan.description = [ ...
    'Nested study-stratified random subsamples. Within each repetition, ' ...
    'the smaller sample is contained in every larger sample.'];
plan.standard = make_cohort_plan( ...
    'BDE/GDE/OBN/AED/VRS', standardSampleSizes, ...
    standardLabels, 810000, 1:67, plan.no_repetitions);
plan.meq30 = make_cohort_plan( ...
    'MEQ30 rows 20:67', meqSampleSizes, ...
    meqLabels, 820000, 20:67, plan.no_repetitions);

verify_plan(plan, standardLabels, meqLabels, standardSampleSizes, ...
    meqSampleSizes, noRepetitions, rngSettings);
save(planPath, 'plan', '-v7');
end


function cohort = make_cohort_plan(name, sampleSizes, labels, seedBase, ...
        originalRows, noRepetitions)
studyIds = unique(labels)';
cohort.name = name;
cohort.sample_sizes = sampleSizes;
cohort.original_rows = originalRows;
cohort.study_ids = studyIds;
cohort.full_study_counts = counts_by_study(labels);
cohort.retained_study_counts = zeros(numel(sampleSizes), 3);
cohort.sample_seeds = seedBase + (1:noRepetitions);
cohort.sample_indices = cell(numel(sampleSizes), noRepetitions);

for sizeIndex = 1:numel(sampleSizes)
    cohort.retained_study_counts(sizeIndex, :) = ...
        proportional_quotas(sampleSizes(sizeIndex), labels);
end

for repetition = 1:noRepetitions
    rng(cohort.sample_seeds(repetition), 'twister');
    orderedByStudy = cell(1, numel(studyIds));
    for groupIndex = 1:numel(studyIds)
        candidates = find(labels == studyIds(groupIndex));
        orderedByStudy{groupIndex} = ...
            candidates(randperm(numel(candidates)));
    end

    for sizeIndex = 1:numel(sampleSizes)
        quotas = cohort.retained_study_counts(sizeIndex, :);
        selected = [];
        for groupIndex = 1:numel(studyIds)
            studyId = studyIds(groupIndex);
            takeN = quotas(studyId);
            selected = [selected; ...
                orderedByStudy{groupIndex}(1:takeN)]; %#ok<AGROW>
        end
        cohort.sample_indices{sizeIndex, repetition} = sort(selected(:));
    end
end
end


function quotas = proportional_quotas(sampleN, labels)
fullCounts = counts_by_study(labels);
expected = sampleN * fullCounts / numel(labels);
quotas = floor(expected);
remaining = sampleN - sum(quotas);
if remaining > 0
    fractions = expected - quotas;
    ranking = sortrows([-fractions(:), (1:3)'], [1 2]);
    for index = 1:remaining
        quotas(ranking(index, 2)) = quotas(ranking(index, 2)) + 1;
    end
end
assert(all(quotas <= fullCounts) && sum(quotas) == sampleN, ...
    'CPM:LearningCurve:QuotaFailure', ...
    'Could not construct proportional quotas for N=%d.', sampleN);
end


function counts = counts_by_study(labels)
counts = arrayfun(@(studyId) sum(labels == studyId), 1:3);
end


function verify_plan(plan, standardLabels, meqLabels, standardSampleSizes, ...
        meqSampleSizes, noRepetitions, rngSettings)
expectedSeeds = rngSettings.cv_seed + (0:(noRepetitions - 1));
assert(plan.version == 3 && plan.no_repetitions == noRepetitions && ...
    isequal(plan.cv_seeds, expectedSeeds) && ...
    strcmp(plan.cv_generator, rngSettings.cv_generator), ...
    'CPM:LearningCurve:PlanVersion', ...
    'The saved learning-curve plan has unexpected global settings.');
verify_cohort(plan.standard, standardLabels, standardSampleSizes, ...
    1:67, plan.no_repetitions);
verify_cohort(plan.meq30, meqLabels, meqSampleSizes, ...
    20:67, plan.no_repetitions);
end


function verify_cohort(cohort, labels, expectedSizes, expectedRows, ...
        noRepetitions)
assert(isequal(cohort.sample_sizes, expectedSizes) && ...
    isequal(cohort.original_rows, expectedRows) && ...
    isequal(size(cohort.sample_indices), ...
        [numel(expectedSizes), noRepetitions]), ...
    'CPM:LearningCurve:CohortPlan', ...
    'A saved cohort plan has unexpected dimensions or sample sizes.');
assert(isequal(cohort.full_study_counts, counts_by_study(labels)), ...
    'CPM:LearningCurve:CohortCounts', ...
    'A saved cohort plan does not match the current study labels.');

for repetition = 1:noRepetitions
    previous = [];
    for sizeIndex = 1:numel(expectedSizes)
        selected = cohort.sample_indices{sizeIndex, repetition};
        assert(numel(selected) == expectedSizes(sizeIndex) && ...
            numel(unique(selected)) == numel(selected) && ...
            all(selected >= 1 & selected <= numel(labels)), ...
            'CPM:LearningCurve:InvalidSelection', ...
            'Invalid selection at N=%d, repetition=%d.', ...
            expectedSizes(sizeIndex), repetition);
        assert(all(ismember(previous, selected)), ...
            'CPM:LearningCurve:NonNestedSelection', ...
            'Samples are not nested at N=%d, repetition=%d.', ...
            expectedSizes(sizeIndex), repetition);
        actualCounts = counts_by_study(labels(selected));
        assert(isequal(actualCounts, ...
            cohort.retained_study_counts(sizeIndex, :)), ...
            'CPM:LearningCurve:StudyQuota', ...
            'Study quotas do not match at N=%d, repetition=%d.', ...
            expectedSizes(sizeIndex), repetition);
        previous = selected;
    end
end
end
