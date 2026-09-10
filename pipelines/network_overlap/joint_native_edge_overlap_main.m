function joint_native_edge_overlap_main( ...
        cfg, outputPath, noIterations, maxWorkers)
%JOINT_NATIVE_EDGE_OVERLAP_MAIN Synchronized native-mask CPM overlap test.
%
% The empirical null jointly permutes Freedman-Lane residuals within the
% two outcome-availability blocks, rebuilds ordinary positive >=8/10 CPM
% consensus masks, and compares a mask-size-standardized overlap score.

if nargin < 2 || isempty(outputPath)
    outputPath = fullfile(cfg.output_root, 'network_overlap');
end
if nargin < 3 || isempty(noIterations)
    noIterations = cfg.network_overlap.iterations;
end
if nargin < 4 || isempty(maxWorkers)
    maxWorkers = cfg.runtime.workers;
end
validateattributes(noIterations, {'numeric'}, {'scalar','integer','positive'});
validateattributes(maxWorkers, {'numeric'}, {'scalar','integer','>=',0});
if ~isfolder(outputPath)
    mkdir(outputPath);
end

pipelineDir = fileparts(mfilename('fullpath'));
repositoryRoot = fileparts(fileparts(pipelineDir));
inputPath = cfg.input_root;
connectomePath = fullfile(inputPath, 'connectomes', 'LSD_difference.mat');
covariatePath = cfg.model.covariate_file;
groupingPath = cfg.model.grouping_file;
entryScript = [mfilename('fullpath') '.m'];
addpath(fullfile(repositoryRoot, 'scripts', 'common'));
addpath(fullfile(repositoryRoot, 'scripts', 'utilities'));
addpath(fileparts(entryScript));

NO_FOLDS = cfg.model.folds;
CONSENSUS_FRACTION = cfg.model.consensus_fraction;
THRESHOLDS = cfg.network_overlap.thresholds;
TOTAL_EDGES = 86320;
MEQ30_ROWS = 20:67;
SHARED_ROWS = 20:67;
EXTRA_ROWS = 1:19;
config = build_outcome_config(cfg.network_overlap.outcomes);

checkpointPath = fullfile(outputPath, 'joint_edge_overlap_checkpoint.mat');
resultsMatPath = fullfile(outputPath, 'joint_edge_overlap_results.mat');
resultsCsvPath = fullfile(outputPath, 'joint_edge_overlap_results.csv');
nullCsvPath = fullfile(outputPath, 'joint_edge_overlap_null_statistics.csv');
diagnosticsCsvPath = fullfile(outputPath, 'joint_edge_overlap_diagnostics.csv');
summaryPath = fullfile(outputPath, 'joint_edge_overlap_summary.txt');
logPath = fullfile(outputPath, 'joint_edge_overlap_run_log.txt');

settings = struct;
settings.analysis = 'joint_native_cpm_edge_overlap_freedman_lane';
settings.no_iterations = noIterations;
settings.no_folds = NO_FOLDS;
settings.thresholds = THRESHOLDS;
settings.consensus_fraction = CONSENSUS_FRACTION;
settings.consensus_folds_required = ceil(CONSENSUS_FRACTION * NO_FOLDS);
settings.covariates = 'full-rank indicator coding for represented studies';
settings.fold_rng = sprintf("rng(%d,'%s')", ...
    cfg.rng.cv_seed, cfg.rng.cv_generator);
settings.permutation_rng = sprintf( ...
    "RandStream('%s','Seed',iteration+%d)", ...
    cfg.rng.internal_permutation_generator, ...
    cfg.rng.permutation_seed_offset);
settings.permutation_draw_order = ...
    'randperm(stream,48) then randperm(stream,19)';
settings.permutation_blocks = ...
    'rows 20:67 joint across GDE/MEQ30/VRS; rows 1:19 joint across GDE/VRS';
settings.network_direction = 'positive';
settings.connectome_form = 'LSD-PCB difference';
settings.total_possible_edges = TOTAL_EDGES;
settings.mask_policy = ...
    'native >=8/10 consensus; no ranking, padding, filling, or discarding';
settings.overlap_statistic = ...
    '-log10(one-sided hypergeometric P[X>=x | K1,K2,86320])';
settings.empty_mask_rule = 'inner p=1 and score=0';
settings.empirical_p = '(1 + count(null score >= observed score))/(B+1)';
settings.primary_correction = 'BH across three p<.01 outcome pairs';
settings.permissive_role = 'p<.05 descriptive only';

[inputs, code] = cpm_analysis_dependencies(cfg, {pipelineDir}, {connectomePath});
guardSettings = struct('settings',settings, 'outcomes',{cfg.network_overlap.outcomes});
cpm_checkpoint_guard(outputPath, guardSettings, inputs, code);
diary(logPath);
diaryCleanup = onCleanup(@() diary('off'));
fprintf('Joint native-mask CPM edge-overlap analysis\n');
fprintf('Invocation started: %s\n', timestamp_now_local());
fprintf('Output: %s\n', outputPath);
fprintf('Iterations: %d; workers requested: %d\n', noIterations, maxWorkers);
fprintf('Outcomes: GDE, MEQ30, VRS; thresholds: p<.01 and p<.05\n');
fprintf('Native masks only: no filler edges and no discarded permutations.\n');

runRecord = struct;
runRecord.status = 'RUNNING';
runRecord.invocation_started = timestamp_now_local();
runRecord.analysis = settings.analysis;
runRecord.purpose = ...
    'Synchronized native positive-consensus edge overlap for Table S16';
runRecord.entry_script = entryScript;
runRecord.input_path = inputPath;
runRecord.output_path = outputPath;
runRecord.connectome_input = connectomePath;
runRecord.covariates_input_file = covariatePath;
runRecord.covariate_columns = ...
    'GDE/VRS: study_2, study_3, sex, age, mean_FD; MEQ30: study_3, sex, age, mean_FD';
runRecord.behaviors = {config.label};
runRecord.behavior_inputs = {config.preferred_behavior};
runRecord.subject_subsets = ...
    'GDE/VRS rows 1:67; MEQ30 rows 20:67';
runRecord.difference_behavior_fallback = ...
    'prefer *_difference behavior, otherwise standard behavior';
runRecord.correlation_types = 'GDE Spearman; MEQ30/VRS Pearson';
runRecord.connectome_form = settings.connectome_form;
runRecord.network_direction = settings.network_direction;
runRecord.no_iterations = noIterations;
runRecord.internal_folds = NO_FOLDS;
runRecord.edge_thresholds = THRESHOLDS;
runRecord.consensus_fraction = CONSENSUS_FRACTION;
runRecord.consensus_folds_required = settings.consensus_folds_required;
runRecord.permutation_scheme = ...
    'Freedman-Lane residual permutation with synchronized availability blocks';
runRecord.permutation_blocks = settings.permutation_blocks;
runRecord.random_seeds = [settings.fold_rng '; ' ...
    settings.permutation_rng '; ' settings.permutation_draw_order];
runRecord.mask_policy = settings.mask_policy;
runRecord.overlap_statistic = settings.overlap_statistic;
runRecord.empty_mask_rule = settings.empty_mask_rule;
runRecord.empirical_p = settings.empirical_p;
runRecord.multiple_comparison_correction = settings.primary_correction;
runRecord.parallel_workers_requested = maxWorkers;
runRecord.checkpoint_mat = checkpointPath;
runRecord.results_mat = resultsMatPath;
runRecord.results_csv = resultsCsvPath;
runRecord.null_statistics_csv = nullCsvPath;
runRecord.diagnostics_csv = diagnosticsCsvPath;
runRecord.results_text = summaryPath;
runRecord.canonical_result_files = {resultsMatPath, resultsCsvPath, ...
    nullCsvPath, diagnosticsCsvPath, summaryPath};
runRecord.completed_result_rows = 0;
runRecord.resumed_from_existing_checkpoint = isfile(checkpointPath);
runRecordPath = save_analysis_run_record(outputPath, '', runRecord);

runTimer = tic;
if maxWorkers > 0
    pool = gcp('nocreate');
    if isempty(pool)
        parpool('Processes', maxWorkers);
        pool = gcp('nocreate');
    else
        assert(pool.NumWorkers == maxWorkers, ...
            'CPM:JointOverlap:UnexpectedPoolSize', ...
            'Existing pool has %d workers; requested %d.', ...
            pool.NumWorkers, maxWorkers);
    end
    activeWorkers = pool.NumWorkers;
else
    activeWorkers = 0;
end
runRecord.parallel_workers_active = activeWorkers;
runRecordPath = save_analysis_run_record(outputPath, runRecordPath, runRecord);

matsFull = load_first_var_local(connectomePath);
covarsFull = load_first_var_local(covariatePath);
assert(isequal(size(matsFull), [416,416,67]) && all(isfinite(matsFull(:))), ...
    'CPM:JointOverlap:InvalidConnectomes', ...
    'LSD_difference must be finite 416-by-416-by-67.');
assert(isequal(size(covarsFull), [67,5]) && all(isfinite(covarsFull(:))), ...
    'CPM:JointOverlap:InvalidCovariates', ...
    'covars_indicator.mat must be finite 67-by-5.');
groupingFull = load_first_var_local(groupingPath);
assert(isequal(size(groupingFull), [67,4]), ...
    'CPM:JointOverlap:InvalidGrouping', 'covars.mat must be 67-by-4.');
studyLabelsFull = groupingFull(:,1);

[behaviors, behaviorDetails] = load_behaviors(config, inputPath);
[covariates, covariateDetails] = build_outcome_covariates( ...
    covarsFull, studyLabelsFull, MEQ30_ROWS);
[behaviorColumns, permutationDetails] = build_joint_fl_columns( ...
    behaviors, covariates, covariateDetails, noIterations, ...
    SHARED_ROWS, EXTRA_ROWS, cfg.rng);
verify_joint_permutations(permutationDetails, noIterations);

fprintf('\nReconstructing observed masks and validating fixed counts...\n');
[observed, edgeLayout] = reconstruct_observed_masks( ...
    config, behaviors, matsFull, covariates, ...
    NO_FOLDS, CONSENSUS_FRACTION, THRESHOLDS, MEQ30_ROWS, ...
    maxWorkers, cfg.rng);
assert(edgeLayout.no_edges == TOTAL_EDGES && edgeLayout.uses_unique_edges, ...
    'CPM:JointOverlap:UnexpectedEdgeLayout', ...
    'Expected 86,320 unique undirected edges.');

checkpoint = initialize_or_load_checkpoint( ...
    checkpointPath, settings, observed, behaviorDetails, permutationDetails);
for outcomeIndex = 1:numel(config)
    if ~isempty(checkpoint.outcomes{outcomeIndex})
        fprintf('\nSkipping completed outcome: %s\n', config(outcomeIndex).label);
        continue
    end
    current = config(outcomeIndex);
    if current.use_meq30_subset
        subjectRows = MEQ30_ROWS;
    else
        subjectRows = 1:67;
    end
    fprintf('\nOutcome %d/3: %s (n=%d, %s)\n', outcomeIndex, ...
        current.label, numel(subjectRows), current.correlation_type);
    rng(cfg.rng.cv_seed, cfg.rng.cv_generator);
    fixedPartition = cvpartition(numel(subjectRows), 'KFold', NO_FOLDS);
    raw = cpm_build_native_consensus_sets( ...
        matsFull(:,:,subjectRows), behaviorColumns{outcomeIndex}, ...
        covariates{outcomeIndex}, current.correlation_type, ...
        fixedPartition, THRESHOLDS, CONSENSUS_FRACTION, maxWorkers);

    for thresholdIndex = 1:numel(THRESHOLDS)
        reconstructed = false(TOTAL_EDGES,1);
        reconstructed(double(raw.positive_sets{1,thresholdIndex})) = true;
        assert(isequal(reconstructed, observed(outcomeIndex).edge_masks{thresholdIndex}), ...
            'CPM:JointOverlap:ObservedReconstructionMismatch', ...
            'Joint engine did not reconstruct observed %s at p<%.2f.', ...
            current.label, THRESHOLDS(thresholdIndex));
    end

    outcome = struct;
    outcome.behavior = current.label;
    outcome.subject_rows = subjectRows;
    outcome.no_subjects = numel(subjectRows);
    outcome.correlation_type = current.correlation_type;
    outcome.thresholds = THRESHOLDS;
    outcome.observed_counts = double(raw.positive_counts(1,:));
    outcome.null_counts = double(raw.positive_counts(2:end,:));
    outcome.null_positive_sets = raw.positive_sets(2:end,:);
    outcome.fixed_fold_test_indices = fold_test_indices(fixedPartition);
    outcome.permutation_seeds = (201:(noIterations+200))';
    outcome.mask_policy = raw.mask_policy;
    checkpoint.outcomes{outcomeIndex} = outcome;
    checkpoint.completed_behaviors{end+1} = current.label;
    checkpoint.last_completed_at = timestamp_now_local();
    cpm_atomic_save(checkpointPath, struct('checkpoint',checkpoint));
    runRecord.completed_result_rows = ...
        sum(~cellfun(@isempty, checkpoint.outcomes));
    runRecordPath = save_analysis_run_record( ...
        outputPath, runRecordPath, runRecord);
    fprintf('Saved native-mask checkpoint for %s.\n', current.label);
    clear raw
end

[resultsTable, nullTable, diagnosticsTable] = build_result_tables( ...
    checkpoint, observed, noIterations, TOTAL_EDGES);
writetable(resultsTable, resultsCsvPath);
writetable(nullTable, nullCsvPath);
writetable(diagnosticsTable, diagnosticsCsvPath);
analysisSettings = settings;
outcomeResults = checkpoint.outcomes;
observedMasks = observed;
save(resultsMatPath, 'analysisSettings', 'outcomeResults', ...
    'observedMasks', 'resultsTable', 'nullTable', 'diagnosticsTable', ...
    'behaviorDetails', 'permutationDetails', '-v7.3');
write_summary(summaryPath, resultsTable, diagnosticsTable, settings);

runRecord.status = 'COMPLETE';
runRecord.analysis_finished = timestamp_now_local();
runRecord.runtime_seconds = toc(runTimer);
runRecord.completed_result_rows = height(resultsTable);
runRecord.completed_null_rows = height(nullTable);
runRecordPath = save_analysis_run_record( ...
    outputPath, runRecordPath, runRecord);
save(fullfile(outputPath, 'run_record_state.mat'), ...
    'runRecord', 'runRecordPath');
fprintf('\nAnalysis COMPLETE: %s\n', timestamp_now_local());
fprintf('Runtime: %.3f seconds\n', runRecord.runtime_seconds);
fprintf('Results: %s\n', resultsCsvPath);
clear diaryCleanup
diary off
end


function config = build_outcome_config(requestedOutcomes)
available = struct( ...
    'label', {'GDE','MEQ30','VRS'}, ...
    'preferred_behavior', {'LSD_GDE_difference', ...
        'LSD_MEQ30_difference','LSD_VRS_difference'}, ...
    'fallback_behavior', {'LSD_GDE','LSD_MEQ30','LSD_VRS'}, ...
    'correlation_type', {'Spearman','Pearson','Pearson'}, ...
    'use_meq30_subset', {false,true,false});
requestedOutcomes = cellstr(string(requestedOutcomes));
config = available(ismember({available.label}, requestedOutcomes));
assert(isequal({config.label}, requestedOutcomes) && numel(config) == 3, ...
    'CPM:JointOverlap:ConfiguredOutcomes', ...
    'Network overlap requires GDE, MEQ30, and VRS in that order.');
end


function [behaviors, details] = load_behaviors(config, inputPath)
behaviors = cell(1,numel(config));
details = repmat(struct('label','', 'source','', 'path','', ...
    'fallback_used',false), numel(config),1);
for index = 1:numel(config)
    [behaviors{index}, source, fallbackUsed, pathName] = ...
        load_behavior_with_fallback(inputPath, ...
        config(index).preferred_behavior, config(index).fallback_behavior);
    details(index).label = config(index).label;
    details(index).source = source;
    details(index).path = pathName;
    details(index).fallback_used = fallbackUsed;
end
end


function [covariates, details] = build_outcome_covariates( ...
        covarsFull, studyLabelsFull, meqRows)
rows = {1:67, meqRows, 1:67};
covariates = cell(1, 3);
details = repmat(struct('columns',{{}}, 'reference_study',NaN, ...
    'design_rank',NaN), 1, 3);
for index = 1:3
    currentRows = rows{index};
    [covariates{index}, names, referenceStudy] = ...
        cpm_select_covariates(covarsFull(currentRows, :), ...
            studyLabelsFull(currentRows));
    details(index).columns = names;
    details(index).reference_study = referenceStudy;
    details(index).design_rank = rank( ...
        [ones(numel(currentRows), 1), covariates{index}]);
end
end


function [columns, details] = build_joint_fl_columns( ...
        behaviors, covariates, covariateDetails, noIterations, ...
        sharedRows, extraRows, rngSettings)
fits = cell(1,3);
residuals = cell(1,3);
rows = {1:67, 20:67, 1:67};
columns = cell(1,3);
for index = 1:3
    currentRows = rows{index};
    outcome = behaviors{index};
    design = [ones(numel(currentRows),1), covariates{index}];
    coefficient = lsqminnorm(design, outcome);
    fits{index} = design * coefficient;
    residuals{index} = outcome - fits{index};
    columns{index} = zeros(numel(currentRows), noIterations+1, 'like', outcome);
    columns{index}(:,1) = outcome;
end

sharedOrders = zeros(noIterations, numel(sharedRows), 'uint8');
extraOrders = zeros(noIterations, numel(extraRows), 'uint8');
for iteration = 1:noIterations
    stream = RandStream(rngSettings.internal_permutation_generator, ...
        'Seed', iteration + rngSettings.permutation_seed_offset);
    pi48 = randperm(stream, numel(sharedRows));
    pi19 = randperm(stream, numel(extraRows));
    order67 = [pi19, numel(extraRows) + pi48];
    columns{1}(:,iteration+1) = fits{1} + residuals{1}(order67);
    columns{2}(:,iteration+1) = fits{2} + residuals{2}(pi48);
    columns{3}(:,iteration+1) = fits{3} + residuals{3}(order67);
    sharedOrders(iteration,:) = uint8(pi48);
    extraOrders(iteration,:) = uint8(pi19);
end
details = struct;
details.scale = 'raw outcome';
details.design_columns = {covariateDetails.columns};
details.study_references = [covariateDetails.reference_study];
details.generator = rngSettings.internal_permutation_generator;
details.seed_range = [1,noIterations] + ...
    rngSettings.permutation_seed_offset;
details.shared_block_rows = sharedRows;
details.extra_block_rows = extraRows;
details.shared_orders = sharedOrders;
details.extra_orders = extraOrders;
details.draw_order = 'shared 48 first, extra 19 second';
end


function verify_joint_permutations(details, noIterations)
assert(isequal(size(details.shared_orders), [noIterations,48]) && ...
    isequal(size(details.extra_orders), [noIterations,19]), ...
    'CPM:JointOverlap:PermutationShapeMismatch', ...
    'Stored block permutations have unexpected dimensions.');
for iteration = 1:noIterations
    assert(isequal(sort(double(details.shared_orders(iteration,:))), 1:48), ...
        'CPM:JointOverlap:InvalidSharedPermutation', ...
        'Shared-block permutation %d is not a bijection.', iteration);
    assert(isequal(sort(double(details.extra_orders(iteration,:))), 1:19), ...
        'CPM:JointOverlap:InvalidExtraPermutation', ...
        'Extra-block permutation %d is not a bijection.', iteration);
end
fprintf('Verified %d synchronized two-block permutation mappings.\n', noIterations);
end


function [observed, referenceLayout] = reconstruct_observed_masks( ...
        config, behaviors, matsFull, covariates, ...
        noFolds, consensusFraction, thresholds, meqRows, maxWorkers, ...
        rngSettings)
observed = repmat(struct('behavior','', 'thresholds',thresholds, ...
    'counts',[], 'edge_masks',{{}}), 3, 1);
referenceLayout = [];
for outcomeIndex = 1:3
    if config(outcomeIndex).use_meq30_subset
        subjectRows = meqRows;
    else
        subjectRows = 1:67;
    end
    rng(rngSettings.cv_seed, rngSettings.cv_generator);
    partition = cvpartition(numel(subjectRows), 'KFold', noFolds);
    edgeMasks = cell(1,2);
    counts = zeros(1,2);
    for thresholdIndex = 1:2
        result = cpm_build_consensus_masks( ...
            matsFull(:,:,subjectRows), behaviors{outcomeIndex}, ...
            covariates{outcomeIndex}, config(outcomeIndex).correlation_type, ...
            partition, thresholds(thresholdIndex), consensusFraction, maxWorkers);
        if isempty(referenceLayout)
            referenceLayout = result.layout;
        else
            assert(isequal(result.layout.upper_index, referenceLayout.upper_index), ...
                'CPM:JointOverlap:EdgeLayoutChanged', ...
                'Unique-edge layout changed between outcomes.');
        end
        edgeMasks{thresholdIndex} = logical( ...
            result.observed_pos_mask(referenceLayout.upper_index));
        counts(thresholdIndex) = nnz(edgeMasks{thresholdIndex});
    end
    observed(outcomeIndex).behavior = config(outcomeIndex).label;
    observed(outcomeIndex).counts = counts;
    observed(outcomeIndex).edge_masks = edgeMasks;
    fprintf('  %s observed counts: %s\n', ...
        config(outcomeIndex).label, mat2str(counts));
end
end


function checkpoint = initialize_or_load_checkpoint( ...
        pathName, settings, observed, behaviorDetails, permutationDetails)
if ~isfile(pathName)
    checkpoint = struct;
    checkpoint.settings = settings;
    checkpoint.observed = observed;
    checkpoint.behavior_details = behaviorDetails;
    checkpoint.permutation_details = permutationDetails;
    checkpoint.outcomes = cell(1,3);
    checkpoint.completed_behaviors = {};
    checkpoint.created_at = timestamp_now_local();
    checkpoint.last_completed_at = '';
    cpm_atomic_save(pathName, struct('checkpoint',checkpoint));
    return
end
loaded = load(pathName, 'checkpoint');
checkpoint = loaded.checkpoint;
assert(isequaln(checkpoint.settings, settings), ...
    'CPM:JointOverlap:CheckpointSettingsMismatch', ...
    'Checkpoint settings differ; use a clean output directory.');
for index = 1:3
    assert(isequal(checkpoint.observed(index).edge_masks, ...
        observed(index).edge_masks), ...
        'CPM:JointOverlap:CheckpointObservedMismatch', ...
        'Observed masks differ from the checkpoint.');
end
end


function [results, nullTable, diagnostics] = build_result_tables( ...
        checkpoint, observed, noIterations, totalEdges)
pairs = [1,2;1,3;2,3];
pairNames = ["GDE_MeQ30";"GDE_VRS";"MEQ30_VRS"];
thresholds = checkpoint.settings.thresholds;
nullTable = table((1:noIterations)', 'VariableNames', {'Permutation'});
rows = cell(6,16);
diagRows = cell(6,11);
rowIndex = 0;
for thresholdIndex = 1:2
    for pairIndex = 1:3
        a = pairs(pairIndex,1);
        b = pairs(pairIndex,2);
        maskA = observed(a).edge_masks{thresholdIndex};
        maskB = observed(b).edge_masks{thresholdIndex};
        sizeA = nnz(maskA);
        sizeB = nnz(maskB);
        shared = nnz(maskA & maskB);
        [innerP, observedScore] = overlap_score( ...
            sizeA, sizeB, shared, totalEdges);
        expected = sizeA * sizeB / totalEdges;
        ratio = shared / expected;
        dice = 2 * shared / (sizeA + sizeB);

        nullSizeA = checkpoint.outcomes{a}.null_counts(:,thresholdIndex);
        nullSizeB = checkpoint.outcomes{b}.null_counts(:,thresholdIndex);
        nullShared = zeros(noIterations,1);
        nullInnerP = ones(noIterations,1);
        nullScore = zeros(noIterations,1);
        for iteration = 1:noIterations
            setA = checkpoint.outcomes{a}.null_positive_sets{iteration,thresholdIndex};
            setB = checkpoint.outcomes{b}.null_positive_sets{iteration,thresholdIndex};
            nullShared(iteration) = numel(intersect(setA,setB));
            [nullInnerP(iteration), nullScore(iteration)] = overlap_score( ...
                nullSizeA(iteration), nullSizeB(iteration), ...
                nullShared(iteration), totalEdges);
        end
        permutationP = (1 + sum(nullScore >= observedScore)) / ...
            (noIterations + 1);
        role = ternary(thresholdIndex == 1, "primary", "descriptive");
        reportedP = ternary(thresholdIndex == 1, permutationP, NaN);

        rowIndex = rowIndex + 1;
        rows(rowIndex,:) = {thresholds(thresholdIndex), ...
            string(observed(a).behavior), string(observed(b).behavior), ...
            sizeA, sizeB, shared, expected, ratio, dice, innerP, ...
            observedScore, reportedP, NaN, ...
            "joint native-mask Freedman-Lane", role, totalEdges};

        suffix = strrep(sprintf('p%.2f',thresholds(thresholdIndex)),'.','');
        stem = sprintf('%s_%s', pairNames(pairIndex), suffix);
        nullTable.([stem '_MaskSizeA']) = nullSizeA;
        nullTable.([stem '_MaskSizeB']) = nullSizeB;
        nullTable.([stem '_SharedEdges']) = nullShared;
        nullTable.([stem '_InnerP']) = nullInnerP;
        nullTable.([stem '_Score']) = nullScore;
        nullTable.([stem '_EitherMaskEmpty']) = ...
            (nullSizeA == 0 | nullSizeB == 0);

        sizeAQuantiles = prctile(nullSizeA,[2.5,50,97.5]);
        sizeBQuantiles = prctile(nullSizeB,[2.5,50,97.5]);
        scoreQuantiles = prctile(nullScore,[2.5,50,97.5]);
        diagRows(rowIndex,:) = {thresholds(thresholdIndex), ...
            string(observed(a).behavior), string(observed(b).behavior), ...
            sum(nullSizeA == 0), sum(nullSizeB == 0), ...
            sizeAQuantiles(1), sizeAQuantiles(2), sizeAQuantiles(3), ...
            sizeBQuantiles(1), sizeBQuantiles(2), sizeBQuantiles(3)};
        fprintf(['p<%.2f %s-%s: score %.6f, exceedances %d, p %.6f; ' ...
            'null score [2.5/50/97.5] %.3f/%.3f/%.3f\n'], ...
            thresholds(thresholdIndex), observed(a).behavior, ...
            observed(b).behavior, observedScore, ...
            sum(nullScore >= observedScore), permutationP, scoreQuantiles);
    end
end
results = cell2table(rows, 'VariableNames', { ...
    'EdgeThreshold','OutcomeA','OutcomeB','MaskSizeA','MaskSizeB', ...
    'SharedEdges','ExpectedSharedEdges','ObservedExpectedRatio','Dice', ...
    'ObservedInnerHypergeometricP','ObservedScore','PermutationP', ...
    'BH_Q_Primary','NullMethod','AnalysisRole','TotalPossibleEdges'});
primary = results.EdgeThreshold == 0.01;
results.BH_Q_Primary(primary) = bh_adjust(results.PermutationP(primary));
diagnostics = cell2table(diagRows, 'VariableNames', { ...
    'EdgeThreshold','OutcomeA','OutcomeB','ZeroMasksA','ZeroMasksB', ...
    'MaskSizeA_P2_5','MaskSizeA_Median','MaskSizeA_P97_5', ...
    'MaskSizeB_P2_5','MaskSizeB_Median','MaskSizeB_P97_5'});
end


function [pValue, score] = overlap_score(sizeA, sizeB, shared, totalEdges)
if sizeA == 0 || sizeB == 0 || shared == 0
    pValue = 1;
    score = 0;
    return
end
pValue = hygecdf(shared-1, totalEdges, sizeA, sizeB, 'upper');
pValue = max(pValue, realmin('double'));
score = -log10(pValue);
assert(isfinite(score), 'CPM:JointOverlap:NonFiniteScore', ...
    'Overlap score is not finite.');
end


function adjusted = bh_adjust(pValues)
pValues = pValues(:);
[sortedP, order] = sort(pValues, 'ascend');
n = numel(pValues);
sortedQ = sortedP .* n ./ (1:n)';
sortedQ = flipud(cummin(flipud(sortedQ)));
sortedQ = min(sortedQ, 1);
adjusted = zeros(size(pValues));
adjusted(order) = sortedQ;
end


function write_summary(pathName, results, diagnostics, settings)
fileId = fopen(pathName, 'wt');
assert(fileId ~= -1, 'CPM:JointOverlap:CannotWriteSummary', ...
    'Could not write %s.', pathName);
cleanup = onCleanup(@() fclose(fileId));
fprintf(fileId, 'Joint native-mask CPM edge-overlap analysis\n');
fprintf(fileId, 'Completed: %s\n\n', timestamp_now_local());
fprintf(fileId, ['Positive, covariate-adjusted LSD-PCB difference masks; ' ...
    '10 fixed folds; native >=8/10 consensus.\n']);
fprintf(fileId, 'Synchronized Freedman-Lane permutations: %d.\n', ...
    settings.no_iterations);
fprintf(fileId, ['Shared rows 20:67 permuted jointly across all outcomes; ' ...
    'rows 1:19 jointly across GDE/VRS.\n']);
fprintf(fileId, ['Statistic: -log10 of the mask-specific one-sided ' ...
    'hypergeometric tail; final p is empirical.\n']);
fprintf(fileId, ['No mask was resized, padded, filled, ranked, or discarded. ' ...
    'BH correction covers the three p<.01 pairs only.\n\n']);
for row = 1:height(results)
    if results.EdgeThreshold(row) == 0.01
        inference = sprintf(', permutation p=%.6f, q=%.6f', ...
            results.PermutationP(row), results.BH_Q_Primary(row));
    else
        inference = ' (descriptive; no p/q reported)';
    end
    fprintf(fileId, ['p<%.2f %s-%s: masks %d/%d, shared %d, expected %.3f, ' ...
        'ratio %.3f, Dice %.4f, inner p %.6g, score %.6f%s\n'], ...
        results.EdgeThreshold(row), results.OutcomeA(row), ...
        results.OutcomeB(row), results.MaskSizeA(row), ...
        results.MaskSizeB(row), results.SharedEdges(row), ...
        results.ExpectedSharedEdges(row), ...
        results.ObservedExpectedRatio(row), results.Dice(row), ...
        results.ObservedInnerHypergeometricP(row), ...
        results.ObservedScore(row), inference);
end
fprintf(fileId, '\nZero-mask diagnostics by pair/threshold:\n');
for row = 1:height(diagnostics)
    fprintf(fileId, 'p<%.2f %s-%s: A=%d, B=%d.\n', ...
        diagnostics.EdgeThreshold(row), diagnostics.OutcomeA(row), ...
        diagnostics.OutcomeB(row), diagnostics.ZeroMasksA(row), ...
        diagnostics.ZeroMasksB(row));
end
clear cleanup
end


function [behavior, source, fallbackUsed, pathName] = ...
        load_behavior_with_fallback(inputPath, preferred, fallback)
preferredPath = fullfile(inputPath, 'behav', [preferred '.mat']);
fallbackPath = fullfile(inputPath, 'behav', [fallback '.mat']);
if isfile(preferredPath)
    pathName = preferredPath;
    source = preferred;
    fallbackUsed = false;
else
    assert(isfile(fallbackPath), 'CPM:JointOverlap:MissingBehavior', ...
        'Neither %s nor %s exists.', preferredPath, fallbackPath);
    pathName = fallbackPath;
    source = fallback;
    fallbackUsed = true;
end
behavior = load_first_var_local(pathName);
behavior = behavior(:);
assert(all(isfinite(behavior)), 'CPM:JointOverlap:InvalidBehavior', ...
    '%s contains NaN or Inf.', pathName);
end


function value = load_first_var_local(pathName)
loaded = load(pathName);
names = fieldnames(loaded);
assert(~isempty(names), 'CPM:JointOverlap:EmptyMat', ...
    'No variable found in %s.', pathName);
value = loaded.(names{1});
end


function indices = fold_test_indices(partition)
indices = cell(partition.NumTestSets,1);
for fold = 1:partition.NumTestSets
    indices{fold} = find(test(partition,fold));
end
end


function out = ternary(condition, ifTrue, ifFalse)
if condition
    out = ifTrue;
else
    out = ifFalse;
end
end


function value = timestamp_now_local()
value = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS'));
end
