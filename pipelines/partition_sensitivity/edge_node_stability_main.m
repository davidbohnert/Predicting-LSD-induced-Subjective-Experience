function edge_node_stability_main(cfg, cvSeeds, outputPath, maxWorkers)
%EDGE_NODE_STABILITY_MAIN Edge and node selection stability across partitions.
%
% Records, for the paper's primary specification (LSD-placebo difference
% connectome, indicator-coded covariate adjustment, positive network), how
% often every one of the 86,320 undirected edges is selected across 100
% ten-fold cross-validation partitions -- the same partitions (CV seeds
% 123:222) whose prediction distributions appear in Figure S1.
%
% Four quantities are stored per edge:
%   fixedPartitionFoldCount  0-10   folds of the manuscript's seed-123 partition
%   partitionConsensusCount  0-100  partitions meeting the >=8/10-fold rule
%   foldSelectionCount       0-1000 all training folds (100 partitions x 10)
%   meanEdgeR                       mean edge-behavior partial correlation over
%                                   all training folds, whether or not the edge
%                                   was selected in them
%
% Feature selection is identical to the publication pipeline: it calls the
% same CPM_PREPARE_EDGE_LAYOUT, CPM_PREPARE_CORRELATION_CACHE and
% CPM_CACHED_EDGE_SELECTION used by CPM_BUILD_CONSENSUS_MASKS. No
% permutation testing is involved. The seed-123 slice is asserted against
% a direct reconstruction of the fixed-partition mask before saving.
%
% Full run:
%   edge_node_stability_main(cfg, 123:222, outputPath, 5)
%
% Cheap serial smoke test:
%   edge_node_stability_main(cfg, 123:125, outputPath, 0)

if nargin < 2 || isempty(cvSeeds), cvSeeds = cfg.partition.seeds; end

entryScript = [mfilename('fullpath') '.m'];
scriptDir = fileparts(entryScript);
repositoryRoot = fileparts(fileparts(scriptDir));
if nargin < 3 || isempty(outputPath)
    outputPath = fullfile(cfg.output_root, 'partition_sensitivity', 'stability');
end
if nargin < 4 || isempty(maxWorkers), maxWorkers = cfg.runtime.workers; end

validateattributes(cvSeeds, {'numeric'}, ...
    {'vector', 'integer', 'positive'}, mfilename, 'cvSeeds');
validateattributes(maxWorkers, {'numeric'}, ...
    {'scalar', 'integer', '>=', 0}, mfilename, 'maxWorkers');
cvSeeds = double(cvSeeds(:))';
outputPath = char(outputPath);

codebaseScripts = fullfile(repositoryRoot, 'scripts');
addpath(fullfile(codebaseScripts, 'common'));
addpath(fullfile(codebaseScripts, 'utilities'));

inputPath = cfg.input_root;
connectomePath = fullfile(inputPath, 'connectomes', 'LSD_difference.mat');
covariatePath = cfg.model.covariate_file;
atlasPath = cfg.atlas_network_file;
if ~isfolder(outputPath), mkdir(outputPath); end


summaryCsv = fullfile(outputPath, 'edge_node_stability_summary.csv');
summaryTxt = fullfile(outputPath, 'edge_node_stability_summary.txt');
runLog = fullfile(outputPath, 'edge_node_stability_run_log.txt');

NO_FOLDS = cfg.model.folds;
P_THRESHOLD = cfg.model.edge_threshold;
CONSENSUS_FRACTION = cfg.model.consensus_fraction;
CONSENSUS_FOLDS = ceil(CONSENSUS_FRACTION * NO_FOLDS);
MEQ30_ROWS = 20:67;
NO_NODES = 416;
TOTAL_POSSIBLE_EDGES = NO_NODES * (NO_NODES - 1) / 2;
noPartitions = numel(cvSeeds);
noTrainingFolds = noPartitions * NO_FOLDS;
outcomeConfig = build_outcome_config();

settings = struct;
settings.version = 2;
settings.analysis = 'edge_and_node_selection_stability_across_cv_partitions';
settings.outcomes = {outcomeConfig.label};
settings.connectome = 'LSD_difference';
settings.connectome_sha256 = file_sha256(connectomePath);
settings.behavior_requested = {outcomeConfig.preferred_behavior};
settings.behavior_fallback = {outcomeConfig.fallback_behavior};
settings.covariates = ...
    'full-rank indicator coding for the represented studies';
settings.covariate_sha256 = file_sha256(covariatePath);
settings.cv_seeds = cvSeeds;
settings.no_partitions = noPartitions;
settings.no_folds = NO_FOLDS;
settings.no_training_folds = noTrainingFolds;
settings.edge_p_threshold = P_THRESHOLD;
settings.consensus_fraction = CONSENSUS_FRACTION;
settings.consensus_folds_required = CONSENSUS_FOLDS;
settings.correlation_types = {outcomeConfig.correlation_type};
settings.meq30_rows = MEQ30_ROWS;
settings.network_direction = 'positive';
settings.fold_rng = sprintf( ...
    'rng(seed,''%s'') before each cvpartition, seeds %d:%d', ...
    cfg.rng.cv_generator, cvSeeds(1), cvSeeds(end));
settings.permutation_scheme = ...
    'none; observed behavior only, no permutation inference';
settings.network_atlas_file = atlasPath;
settings.network_atlas_sha256 = file_sha256(atlasPath);
settings.total_possible_edges = TOTAL_POSSIBLE_EDGES;
settings.stability_definitions = [ ...
    'fixedPartitionFoldCount: folds of the seed-123 partition selecting ' ...
    'the edge (0-10). partitionConsensusCount: partitions in which the ' ...
    'edge met the >=8/10-fold consensus rule (0-100). ' ...
    'foldSelectionCount: all training folds selecting the edge ' ...
    '(0-100 partitions x 10 folds). meanEdgeR: mean edge-behavior ' ...
    'partial correlation across all training folds, using the ' ...
    'outcome''s own correlation method.'];

[inputs, code] = cpm_analysis_dependencies(cfg, {scriptDir}, ...
    {connectomePath, atlasPath});
guardSettings = struct('seeds',cvSeeds, 'folds',NO_FOLDS, ...
    'threshold',P_THRESHOLD,'consensus',CONSENSUS_FRACTION, ...
    'rng',cfg.rng,'outcomes',{settings.outcomes});
cpm_checkpoint_guard(outputPath, guardSettings, inputs, code);
copyfile(cfg.atlas_parcel_metadata_file, fullfile(outputPath,'atlas_parcel_metadata.csv'));
runTimer = tic;
diary(runLog);
diary on;
diaryCleanup = onCleanup(@() diary('off'));

fprintf('CPM edge and node selection stability across CV partitions\n');
fprintf('Invocation started: %s\n', timestamp_now_local());
fprintf('Output: %s\n', outputPath);
fprintf('Outcomes: %s\n', strjoin({outcomeConfig.label}, ', '));
fprintf('Partitions: %d (CV seeds %d:%d)\n', ...
    noPartitions, cvSeeds(1), cvSeeds(end));
fprintf('Feature selection: p < %.3g; consensus >= %d/%d folds\n', ...
    P_THRESHOLD, CONSENSUS_FOLDS, NO_FOLDS);
fprintf('Workers requested: %d\n', maxWorkers);

runRecord = struct;
runRecord.status = 'RUNNING';
runRecord.invocation_started = timestamp_now_local();
runRecord.analysis = settings.analysis;
runRecord.purpose = ['Supplementary Table S17: most consistently selected ' ...
    'positive edges and candidate high-degree nodes across 100 ten-fold ' ...
    'cross-validation partitions'];
runRecord.entry_script = entryScript;
runRecord.input_path = inputPath;
runRecord.output_path = outputPath;
runRecord.connectome_input = connectomePath;
runRecord.covariates_input_file = covariatePath;
runRecord.covariate_columns = ...
    'GDE/VRS: study_2, study_3, sex, age, mean_FD; MEQ30: study_3, sex, age, mean_FD';
runRecord.behaviors = {outcomeConfig.label};
runRecord.behavior_inputs = {outcomeConfig.preferred_behavior};
runRecord.connectome_form = 'LSD-PCB difference';
runRecord.network_direction = 'positive';
runRecord.no_iterations = 0;
runRecord.no_partitions = noPartitions;
runRecord.internal_folds = NO_FOLDS;
runRecord.edge_threshold = P_THRESHOLD;
runRecord.consensus_fraction = CONSENSUS_FRACTION;
runRecord.consensus_folds_required = CONSENSUS_FOLDS;
runRecord.correlation_types = 'GDE Spearman; MEQ30 Pearson; VRS Pearson';
runRecord.subject_subsets = 'GDE/VRS rows 1:67; MEQ30 rows 20:67';
runRecord.difference_behavior_fallback = ...
    'prefer *_difference behavior, otherwise standard behavior';
runRecord.permutation_scheme = settings.permutation_scheme;
runRecord.random_seeds = settings.fold_rng;
runRecord.network_atlas_file = atlasPath;
runRecord.stability_definitions = settings.stability_definitions;
runRecord.reproduction_check = ...
    'seed-123 mask is rebuilt directly by the publication implementation';
runRecord.parallel_workers_requested = maxWorkers;
runRecord.results_text = summaryTxt;
runRecord.summary_csv = summaryCsv;
runRecord.frequency_mat_pattern = fullfile(outputPath, ...
    'edge_selection_frequencies_<outcome>.mat');
runRecord.canonical_result_files = [{summaryCsv, summaryTxt}, ...
    arrayfun(@(c) frequency_path(outputPath, c.label), outcomeConfig, ...
    'UniformOutput', false)];
runRecord.completed_result_rows = 0;
runRecordPath = save_analysis_run_record(outputPath, '', runRecord);

if maxWorkers > 0
    pool = gcp('nocreate');
    if ~isempty(pool)
        assert(pool.NumWorkers == maxWorkers, ...
            'CPM:EdgeNodeStability:UnexpectedExistingPool', ...
            ['An existing pool has %d workers; requested %d. Stop the ' ...
             'other analysis/pool before this run.'], ...
            pool.NumWorkers, maxWorkers);
    else
        parpool('Processes', maxWorkers);
        pool = gcp('nocreate');
    end
    activeWorkers = pool.NumWorkers;
else
    activeWorkers = 0;
end
runRecord.parallel_workers_active = activeWorkers;
runRecordPath = save_analysis_run_record(outputPath, runRecordPath, runRecord);
fprintf('Parallel workers active: %d\n', activeWorkers);

matsFull = load_first_var_local(connectomePath);
validateattributes(matsFull, {'numeric'}, ...
    {'3d', 'real', 'finite', 'nonempty'});
assert(isequal(size(matsFull), [NO_NODES, NO_NODES, 67]), ...
    'CPM:EdgeNodeStability:UnexpectedConnectomeSize', ...
    'LSD_difference must be 416-by-416-by-67.');
covarsFull = load_first_var_local(covariatePath);
assert(isnumeric(covarsFull) && isequal(size(covarsFull), [67, 5]) && ...
    all(isfinite(covarsFull(:))), ...
    'CPM:EdgeNodeStability:InvalidCovariates', ...
    ['covars_indicator.mat must be finite 67-by-5: ' ...
     '[study_2, study_3, sex, age, mean_FD].']);
groupingFull = load_first_var_local(cfg.model.grouping_file);
assert(isequal(size(groupingFull), [67, 4]), ...
    'CPM:EdgeNodeStability:InvalidGrouping', ...
    'covars.mat must be 67-by-4.');
studyLabelsFull = groupingFull(:, 1);

summaryRows = cell(numel(outcomeConfig), 1);
for outcomeIndex = 1:numel(outcomeConfig)
    config = outcomeConfig(outcomeIndex);

    [behavior, behaviorSource, fallbackUsed, behaviorPath] = ...
        load_behavior_with_fallback(inputPath, ...
        config.preferred_behavior, config.fallback_behavior);
    if config.use_meq30_subset
        subjectRows = MEQ30_ROWS;
    else
        subjectRows = 1:67;
    end
    mats = matsFull(:, :, subjectRows);
    [covars, covariateNames, referenceStudy] = ...
        cpm_select_covariates(covarsFull(subjectRows, :), ...
            studyLabelsFull(subjectRows));
    assert(numel(behavior) == numel(subjectRows), ...
        'CPM:EdgeNodeStability:BehaviorCountMismatch', ...
        '%s contains %d values; expected %d.', ...
        behaviorSource, numel(behavior), numel(subjectRows));

    fprintf('\n========================================\n');
    fprintf('Outcome: %s (n=%d, %s)\n', ...
        config.label, numel(subjectRows), config.correlation_type);
    fprintf('Behavior source: %s\n', behaviorSource);
    fprintf('========================================\n');
    outcomeTimer = tic;

    frequencyPath = frequency_path(outputPath, config.label);
    if isfile(frequencyPath)
        saved = load(frequencyPath, 'foldSelectionCount', ...
            'partitionConsensusCount', 'fixedPartitionFoldCount', 'metadata');
        assert(isequal(saved.metadata.cv_seeds, cvSeeds) && ...
            isequal(saved.metadata.covariate_columns, covariateNames), ...
            'CPM:EdgeNodeStability:CheckpointMismatch', ...
            'Existing %s output uses different settings.', config.label);
        fixedEdgeCount = nnz(saved.fixedPartitionFoldCount >= CONSENSUS_FOLDS);
        summaryRows{outcomeIndex} = make_summary_row(config, noSubjects, ...
            behaviorSource, fixedEdgeCount, saved.foldSelectionCount, ...
            saved.partitionConsensusCount, noPartitions, noTrainingFolds, 0);
        fprintf('Using completed checkpoint: %s\n', frequencyPath);
        runRecord.completed_result_rows = outcomeIndex;
        runRecordPath = save_analysis_run_record( ...
            outputPath, runRecordPath, runRecord);
        continue
    end

    layout = cpm_prepare_edge_layout(mats, config.correlation_type);
    assert(layout.uses_unique_edges && ...
        layout.no_edges == TOTAL_POSSIBLE_EDGES, ...
        'CPM:EdgeNodeStability:UnexpectedEdgeLayout', ...
        'Expected 86,320 unique undirected edges.');

    correlationEdges = layout.correlation_edges;
    noSubjects = numel(subjectRows);
    noEdges = layout.no_edges;

    % One column per partition keeps the PARFOR reduction trivially sliced.
    % uint8 is enough: within one partition an edge is selected 0-10 times.
    partitionFoldCounts = zeros(noEdges, noPartitions, 'uint8');
    partitionEdgeRSums = zeros(noEdges, noPartitions);
    if maxWorkers > 0
        parfor (partitionIndex = 1:noPartitions, maxWorkers)
            [partitionFoldCounts(:, partitionIndex), ...
                partitionEdgeRSums(:, partitionIndex)] = ...
                partition_fold_counts(cvSeeds(partitionIndex), ...
                NO_FOLDS, noSubjects, correlationEdges, behavior, ...
                covars, config.correlation_type, P_THRESHOLD, ...
                cfg.rng.cv_generator); %#ok<PFBNS>
        end
    else
        for partitionIndex = 1:noPartitions
            [partitionFoldCounts(:, partitionIndex), ...
                partitionEdgeRSums(:, partitionIndex)] = ...
                partition_fold_counts(cvSeeds(partitionIndex), ...
                NO_FOLDS, noSubjects, correlationEdges, behavior, ...
                covars, config.correlation_type, P_THRESHOLD, ...
                cfg.rng.cv_generator);
            fprintf('  ...partition %d / %d complete\n', ...
                partitionIndex, noPartitions);
        end
    end

    foldSelectionCount = uint16(sum(double(partitionFoldCounts), 2));
    meanEdgeR = sum(partitionEdgeRSums, 2) / noTrainingFolds;
    assert(all(isfinite(meanEdgeR)) && max(abs(meanEdgeR)) <= 1, ...
        'CPM:EdgeNodeStability:InvalidMeanCorrelation', ...
        'Mean edge correlations must be finite and within [-1, 1].');
    partitionConsensusCount = ...
        uint8(sum(partitionFoldCounts >= CONSENSUS_FOLDS, 2));
    fixedIndex = find(cvSeeds == 123, 1);
    assert(~isempty(fixedIndex), ...
        'CPM:EdgeNodeStability:MissingFixedPartition', ...
        'CV seed 123, the manuscript partition, must be included.');
    fixedPartitionFoldCount = partitionFoldCounts(:, fixedIndex);

    assert(max(foldSelectionCount) <= noTrainingFolds && ...
        max(partitionConsensusCount) <= noPartitions && ...
        max(fixedPartitionFoldCount) <= NO_FOLDS, ...
        'CPM:EdgeNodeStability:CountOutOfRange', ...
        'A stability count exceeds its maximum possible value.');

    fixedEdgeCount = nnz(fixedPartitionFoldCount >= CONSENSUS_FOLDS);
    fprintf('Seed-123 mask rebuilt: %d positive edges\n', ...
        fixedEdgeCount);

    edgeIndex = int32(layout.upper_index);
    metadata = struct;
    metadata.analysis_id = sprintf('LSD_%s_difference_edge_node_stability', ...
        config.label);
    metadata.scale = config.label;
    metadata.matrix_name = 'LSD_difference';
    metadata.connectome_form = 'difference';
    metadata.behavior_requested = config.preferred_behavior;
    metadata.behavior_source = behaviorSource;
    metadata.behavior_fallback_used = fallbackUsed;
    metadata.behavior_path = behaviorPath;
    metadata.behavior_sha256 = file_sha256(behaviorPath);
    metadata.correlation_type = config.correlation_type;
    metadata.covariate_setting = 'with_covariates';
    metadata.covariate_path = covariatePath;
    metadata.covariate_columns = covariateNames;
    metadata.study_reference = referenceStudy;
    metadata.covariate_design_rank = ...
        rank([ones(noSubjects, 1), covars]);
    metadata.subject_rows = subjectRows;
    metadata.no_subjects = noSubjects;
    metadata.node_count = NO_NODES;
    metadata.no_edges = noEdges;
    metadata.no_folds = NO_FOLDS;
    metadata.cv_seeds = cvSeeds;
    metadata.no_partitions = noPartitions;
    metadata.no_training_folds = noTrainingFolds;
    metadata.edge_p_threshold = P_THRESHOLD;
    metadata.consensus_fraction = CONSENSUS_FRACTION;
    metadata.consensus_folds_required = CONSENSUS_FOLDS;
    metadata.network_direction = 'positive';
    metadata.fixed_partition_seed = 123;
    metadata.fixed_partition_edges = fixedEdgeCount;
    metadata.edge_representation = ...
        'unique upper triangle, column-major, edgeIndex into 416x416';
    metadata.stability_definitions = settings.stability_definitions;
    metadata.mean_edge_r_definition = ['mean edge-behavior ' ...
        'partial correlation over all training folds, whether or not the ' ...
        'edge was selected in them'];
    metadata.rng = settings.fold_rng;
    metadata.entry_script = entryScript;
    metadata.created_at = timestamp_now_local();

    save(frequencyPath, 'foldSelectionCount', 'partitionConsensusCount', ...
        'fixedPartitionFoldCount', 'meanEdgeR', 'edgeIndex', 'metadata', ...
        '-v7');
    fprintf('Saved: %s\n', frequencyPath);

    elapsed = toc(outcomeTimer);
    summaryRows{outcomeIndex} = make_summary_row(config, noSubjects, ...
        behaviorSource, fixedEdgeCount, foldSelectionCount, ...
        partitionConsensusCount, noPartitions, noTrainingFolds, elapsed);
    halfPartitions = ceil(noPartitions / 2);
    fprintf(['Edges ever selected: %d; ever reaching consensus: %d; ' ...
        'in >=%d/%d partitions: %d\n'], nnz(foldSelectionCount > 0), ...
        nnz(partitionConsensusCount > 0), halfPartitions, ...
        noPartitions, nnz(partitionConsensusCount >= halfPartitions));
    fprintf('Elapsed: %.1f seconds\n', elapsed);

    runRecord.completed_result_rows = outcomeIndex;
    runRecordPath = save_analysis_run_record(outputPath, runRecordPath, runRecord);
end

summary = cell2table(vertcat(summaryRows{:}), 'VariableNames', ...
    {'Outcome', 'NSubjects', 'CorrelationType', 'BehaviorSource', ...
     'FixedPartitionConsensusEdges', 'EdgesEverSelected', ...
     'EdgesEverReachingConsensus', 'EdgesConsensusInAtLeastHalfPartitions', ...
     'MaxConsensusRecurrence', 'MaxFoldSelectionCount', ...
     'MeanEdgesSelectedPerFold', 'ElapsedSeconds'});
writetable(summary, summaryCsv);
write_summary_text(summaryTxt, summary, settings);

runtimeSeconds = toc(runTimer);
runRecord.status = 'COMPLETE';
runRecord.analysis_finished = timestamp_now_local();
runRecord.runtime_seconds = runtimeSeconds;
runRecord.completed_result_rows = height(summary);
save_analysis_run_record(outputPath, runRecordPath, runRecord);

fprintf('\nAll outcomes complete in %.1f seconds (%.1f minutes).\n', ...
    runtimeSeconds, runtimeSeconds / 60);
fprintf('Summary: %s\n', summaryCsv);
end


function row = make_summary_row(config, noSubjects, behaviorSource, ...
        fixedEdgeCount, foldSelectionCount, partitionConsensusCount, ...
        noPartitions, noTrainingFolds, elapsed)
row = {string(config.label), noSubjects, string(config.correlation_type), ...
    string(behaviorSource), fixedEdgeCount, nnz(foldSelectionCount > 0), ...
    nnz(partitionConsensusCount > 0), ...
    nnz(partitionConsensusCount >= ceil(noPartitions / 2)), ...
    max(double(partitionConsensusCount)), ...
    max(double(foldSelectionCount)), ...
    sum(double(foldSelectionCount)) / noTrainingFolds, elapsed};
end


function [foldCounts, edgeRSum] = partition_fold_counts(seed, noFolds, ...
        noSubjects, correlationEdges, behavior, covars, corrType, ...
        pThreshold, cvGenerator)
%PARTITION_FOLD_COUNTS Per-edge fold-selection counts and correlation sums.
%
% One partition, built with the pipeline's own RNG convention. The body
% mirrors the private SELECT_FOLD_EDGES of CPM_BUILD_CONSENSUS_MASKS so
% that feature selection is bit-for-bit the published implementation; only
% the bookkeeping differs, because the consensus builder discards the
% per-fold votes this table needs.
rng(seed, cvGenerator);
partition = cvpartition(noSubjects, 'KFold', noFolds);
foldCounts = zeros(size(correlationEdges, 1), 1, 'uint8');
edgeRSum = zeros(size(correlationEdges, 1), 1);
for fold = 1:noFolds
    trainIndex = find(training(partition, fold));
    foldEdges = correlationEdges(:, trainIndex)';
    if isempty(covars)
        foldCovars = [];
    else
        foldCovars = covars(trainIndex, :);
    end
    cache = cpm_prepare_correlation_cache(foldEdges, foldCovars, corrType);
    [edgeR, posMasks, ~] = cpm_cached_edge_selection( ...
        cache, behavior(trainIndex), pThreshold);
    foldCounts = foldCounts + uint8(posMasks{1});
    edgeRSum = edgeRSum + edgeR;
end
end


function config = build_outcome_config()
% Identical to BUILD_OUTCOME_CONFIG in
% network_interaction_freedman_lane_main.m, so the two analyses describe
% exactly the same three models.
config = struct( ...
    'label', {'GDE', 'MEQ30', 'VRS'}, ...
    'preferred_behavior', {'LSD_GDE_difference', ...
        'LSD_MEQ30_difference', 'LSD_VRS_difference'}, ...
    'fallback_behavior', {'LSD_GDE', 'LSD_MEQ30', 'LSD_VRS'}, ...
    'correlation_type', {'Spearman', 'Pearson', 'Pearson'}, ...
    'use_meq30_subset', {false, true, false});
end


function pathName = frequency_path(outputPath, label)
pathName = fullfile(outputPath, ...
    sprintf('edge_selection_frequencies_%s.mat', label));
end


function write_summary_text(pathName, summary, settings)
fileId = fopen(pathName, 'wt');
assert(fileId > 0, 'CPM:EdgeNodeStability:CannotWriteSummary', ...
    'Could not write %s', pathName);
cleanup = onCleanup(@() fclose(fileId));
fprintf(fileId, 'CPM edge and node selection stability across CV partitions\n');
fprintf(fileId, 'Generated: %s\n\n', timestamp_now_local());
fprintf(fileId, 'Connectome: %s (LSD-placebo difference)\n', settings.connectome);
fprintf(fileId, 'Covariates: %s\n', settings.covariates);
fprintf(fileId, 'Network direction: %s\n', settings.network_direction);
fprintf(fileId, 'Partitions: %d (CV seeds %d:%d), %d folds each\n', ...
    settings.no_partitions, settings.cv_seeds(1), settings.cv_seeds(end), ...
    settings.no_folds);
fprintf(fileId, 'Training folds in total: %d\n', settings.no_training_folds);
fprintf(fileId, 'Edge threshold: p < %.3g; consensus >= %d of %d folds\n', ...
    settings.edge_p_threshold, settings.consensus_folds_required, ...
    settings.no_folds);
fprintf(fileId, 'Possible undirected edges: %d\n\n', ...
    settings.total_possible_edges);
fprintf(fileId, '%s\n\n', settings.stability_definitions);
for row = 1:height(summary)
    fprintf(fileId, '%s (n=%d, %s, behavior %s)\n', ...
        summary.Outcome(row), summary.NSubjects(row), ...
        summary.CorrelationType(row), summary.BehaviorSource(row));
    fprintf(fileId, '  Seed-123 consensus edges:            %d\n', ...
        summary.FixedPartitionConsensusEdges(row));
    fprintf(fileId, '  Edges selected in at least one fold: %d\n', ...
        summary.EdgesEverSelected(row));
    fprintf(fileId, '  Edges reaching consensus at all:     %d\n', ...
        summary.EdgesEverReachingConsensus(row));
    fprintf(fileId, '  Edges reaching consensus in >=50%%:    %d\n', ...
        summary.EdgesConsensusInAtLeastHalfPartitions(row));
    fprintf(fileId, '  Highest consensus recurrence:        %d\n', ...
        summary.MaxConsensusRecurrence(row));
    fprintf(fileId, '  Highest fold-selection count:        %d\n', ...
        summary.MaxFoldSelectionCount(row));
    fprintf(fileId, '  Mean edges selected per fold:        %.1f\n\n', ...
        summary.MeanEdgesSelectedPerFold(row));
end
end


function [behavior, sourceName, fallbackUsed, behaviorPath] = ...
        load_behavior_with_fallback(inputPath, preferredName, fallbackName)
preferredPath = fullfile(inputPath, 'behav', [preferredName '.mat']);
if isfile(preferredPath)
    sourceName = preferredName;
    behaviorPath = preferredPath;
    fallbackUsed = false;
else
    behaviorPath = fullfile(inputPath, 'behav', [fallbackName '.mat']);
    assert(isfile(behaviorPath), ...
        'CPM:EdgeNodeStability:MissingBehavior', ...
        'Neither preferred nor fallback behavior exists for %s.', preferredName);
    sourceName = fallbackName;
    fallbackUsed = true;
end
behavior = load_first_var_local(behaviorPath);
validateattributes(behavior, {'numeric'}, ...
    {'vector', 'real', 'finite', 'nonempty'});
behavior = behavior(:);
end


function value = load_first_var_local(pathName)
assert(isfile(pathName), 'CPM:EdgeNodeStability:MissingInput', ...
    'Required input is missing: %s', pathName);
loaded = load(pathName);
names = fieldnames(loaded);
assert(isscalar(names), 'CPM:EdgeNodeStability:AmbiguousInput', ...
    'Expected exactly one variable in %s.', pathName);
value = loaded.(names{1});
end


function value = timestamp_now_local()
value = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS'));
end
