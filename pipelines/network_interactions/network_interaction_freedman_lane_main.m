function network_interaction_freedman_lane_main( ...
        cfg, noIterations, outputPath, maxWorkers)
%NETWORK_INTERACTION_FREEDMAN_LANE_MAIN CPM network-pair enrichment inference.
%
% Rebuilds the covariate-adjusted positive 80%-consensus masks for the
% LSD-placebo difference model, then compares raw network-pair edge counts
% with a Freedman-Lane behavioral-permutation null. The observed / expected
% ratio is retained as the displayed effect size. Each permutation keeps
% the connectomes, fixed folds, nuisance signal, edge selection, and
% consensus rule intact.
%
% Full analysis:
%   network_interaction_freedman_lane_main(cfg, 10000, outputPath, 5)
%
% Cheap serial smoke test:
%   network_interaction_freedman_lane_main(cfg, 3, outputPath, 0)

if nargin < 2 || isempty(noIterations)
    noIterations = cfg.network_interactions.iterations;
end

entryScript = [mfilename('fullpath') '.m'];
scriptDir = fileparts(entryScript);
repositoryRoot = fileparts(fileparts(scriptDir));
if nargin < 3 || isempty(outputPath)
    outputPath = fullfile(cfg.output_root, 'network_interactions');
end
if nargin < 4 || isempty(maxWorkers), maxWorkers = cfg.runtime.workers; end

validateattributes(noIterations, {'numeric'}, ...
    {'scalar', 'integer', '>=', 1}, mfilename, 'noIterations');
validateattributes(maxWorkers, {'numeric'}, ...
    {'scalar', 'integer', '>=', 0}, mfilename, 'maxWorkers');
outputPath = char(outputPath);

codebaseScripts = fullfile(repositoryRoot, 'scripts');
addpath(fullfile(codebaseScripts, 'common'));
addpath(fullfile(codebaseScripts, 'utilities'));

inputPath = cfg.input_root;
connectomePath = fullfile(inputPath, 'connectomes', 'LSD_difference.mat');
covariatePath = cfg.model.covariate_file;
groupingPath = cfg.model.grouping_file;
atlasPath = cfg.atlas_network_file;
if ~isfolder(outputPath), mkdir(outputPath); end

summaryCsv = fullfile(outputPath, 'outcome_summary.csv');
pairCsv = fullfile(outputPath, 'network_pair_results.csv');
nullCountCsv = fullfile(outputPath, 'null_consensus_edge_counts.csv');
resultsMat = fullfile(outputPath, 'network_interaction_results.mat');
checkpointPath = fullfile(outputPath, 'network_interaction_checkpoint.mat');
summaryTxt = fullfile(outputPath, 'network_interaction_summary.txt');
runLog = fullfile(outputPath, 'network_interaction_run_log.txt');

NO_FOLDS = cfg.model.folds;
P_THRESHOLD = cfg.model.edge_threshold;
CONSENSUS_FRACTION = cfg.model.consensus_fraction;
CONSENSUS_FOLDS = ceil(CONSENSUS_FRACTION * NO_FOLDS);
MEQ30_ROWS = 20:67;
TOTAL_POSSIBLE_EDGES = 416 * 415 / 2;
outcomeConfig = build_outcome_config(cfg.network_interactions.outcomes);
networkDefinition = load_network_definition(atlasPath, 416);
assert(networkDefinition.no_pairs == cfg.network_interactions.pair_count, ...
    'CPM:NetworkInteraction:ConfiguredPairCount', ...
    'Configured and atlas-derived network-pair counts disagree.');

settings = struct;
settings.version = 4;
settings.analysis = ...
    'network_interaction_freedman_lane_positive_consensus_count_inference';
settings.outcomes = {outcomeConfig.label};
settings.connectome = 'LSD_difference';
settings.connectome_sha256 = file_sha256(connectomePath);
settings.behavior_requested = {outcomeConfig.preferred_behavior};
settings.behavior_fallback = {outcomeConfig.fallback_behavior};
settings.covariates = ...
    'full-rank indicator coding for the represented studies';
settings.covariate_sha256 = file_sha256(covariatePath);
settings.no_iterations = noIterations;
settings.no_folds = NO_FOLDS;
settings.edge_p_threshold = P_THRESHOLD;
settings.consensus_fraction = CONSENSUS_FRACTION;
settings.consensus_folds_required = CONSENSUS_FOLDS;
settings.correlation_types = {outcomeConfig.correlation_type};
settings.meq30_rows = MEQ30_ROWS;
settings.fold_rng = sprintf('rng(%d,''%s'') before each cvpartition', ...
    cfg.rng.cv_seed, cfg.rng.cv_generator);
settings.permutation_rng = sprintf( ...
    'RandStream(''%s'',''Seed'',iteration+%d), iterations 1:%d', ...
    cfg.rng.internal_permutation_generator, ...
    cfg.rng.permutation_seed_offset, noIterations);
settings.permutation_method = ...
    ['Freedman-Lane on raw outcome: y_null = nuisance fitted values + ' ...
     'permuted reduced-model residuals'];
settings.network_atlas_file = atlasPath;
settings.network_atlas_sha256 = file_sha256(atlasPath);
settings.network_names = cellstr(networkDefinition.names);
settings.no_network_pairs = networkDefinition.no_pairs;
settings.total_possible_edges = TOTAL_POSSIBLE_EDGES;
settings.empirical_test = ...
    ['raw network-pair consensus-edge count; one-sided upper tail, ' ...
     '(1 + number of null counts >= observed count) / (permutations + 1)'];
settings.zero_edge_rule = ...
    ['include all null masks, including zero-edge masks, in raw-count ' ...
     'network-pair and total-mask inference'];
settings.display_effect_size = ...
    ['observed / expected edges: pair edge share over overall edge share; ' ...
     '1 = the share an even spread across all connections would give'];
settings.total_mask_test = ...
    ['separate one-sided upper-tail permutation test of total positive ' ...
     'consensus-edge count'];
settings.multiple_comparison = ...
    ['Benjamini-Hochberg FDR across 36 pairs separately within each ' ...
     'outcome; total-mask BH across GDE, MEQ30, and VRS'];
settings.parallel_workers_requested = maxWorkers;

[inputs, code] = cpm_analysis_dependencies(cfg, {scriptDir}, {connectomePath, atlasPath});
guardSettings = struct('model',rmfield(cfg.model,{'covariate_file','grouping_file'}), ...
    'outcomes',{cfg.network_interactions.outcomes}, 'rng',cfg.rng, 'iterations',noIterations);
cpm_checkpoint_guard(outputPath, guardSettings, inputs, code);
runTimer = tic;
diary(runLog);
diary on;
diaryCleanup = onCleanup(@() diary('off'));
fprintf('CPM network-pair enrichment: Freedman-Lane null\n');
fprintf('Invocation started: %s\n', timestamp_now_local());
fprintf('Output: %s\n', outputPath);
fprintf('Outcomes: %s\n', strjoin({outcomeConfig.label}, ', '));
fprintf('Permutations per outcome: %d\n', noIterations);
fprintf('Workers requested: %d\n', maxWorkers);
fprintf('Feature selection: p < %.3g; consensus >= %d/%d folds\n', ...
    P_THRESHOLD, CONSENSUS_FOLDS, NO_FOLDS);
fprintf('Power status is recorded by the terminal launcher.\n');

runRecord = struct;
runRecord.status = 'RUNNING';
runRecord.invocation_started = timestamp_now_local();
runRecord.analysis = settings.analysis;
runRecord.purpose = ['Network-pair concentration of positive CPM consensus ' ...
    'edges using raw pair counts for NCA-style inference and the ' ...
    'observed / expected ratio for display'];
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
runRecord.no_iterations = noIterations;
runRecord.internal_folds = NO_FOLDS;
runRecord.edge_threshold = P_THRESHOLD;
runRecord.consensus_fraction = CONSENSUS_FRACTION;
runRecord.consensus_folds_required = CONSENSUS_FOLDS;
runRecord.correlation_types = ...
    'GDE Spearman; MEQ30 Pearson; VRS Pearson';
runRecord.subject_subsets = ...
    'GDE/VRS rows 1:67; MEQ30 rows 20:67';
runRecord.difference_behavior_fallback = ...
    'prefer *_difference behavior, otherwise standard behavior';
runRecord.permutation_scheme = settings.permutation_method;
runRecord.random_seeds = [settings.fold_rng '; ' settings.permutation_rng];
runRecord.network_atlas_file = atlasPath;
runRecord.network_pair_count_method = ...
    ['unique undirected positive consensus edges counted once; ' ...
     'within-network possible edges n(n-1)/2'];
runRecord.empirical_test = settings.empirical_test;
runRecord.display_effect_size = settings.display_effect_size;
runRecord.total_mask_test = settings.total_mask_test;
runRecord.multiple_comparison = settings.multiple_comparison;
runRecord.parallel_workers_requested = maxWorkers;
runRecord.python_executable = cfg.runtime.python_executable;
runRecord.results_text = summaryTxt;
runRecord.outcome_summary_csv = summaryCsv;
runRecord.network_pair_results_csv = pairCsv;
runRecord.null_edge_counts_csv = nullCountCsv;
runRecord.results_mat = resultsMat;
runRecord.checkpoint_mat = checkpointPath;
runRecord.results_workbook = fullfile(outputPath, ...
    'network_interaction_results.xlsx');
runRecord.figures_directory = fullfile(outputPath, 'figures');
runRecord.canonical_result_files = ...
    {summaryCsv, pairCsv, nullCountCsv, resultsMat, summaryTxt, ...
     runRecord.results_workbook};
runRecord.canonical_figure_files = { ...
    fullfile(outputPath, 'figures', 'network_enrichment_GDE.png'), ...
    fullfile(outputPath, 'figures', 'network_enrichment_MEQ30.png'), ...
    fullfile(outputPath, 'figures', 'network_enrichment_VRS.png')};
runRecord.completed_result_rows = 0;
% One output directory keeps one run record. Re-running to rebuild tables from a
% finished checkpoint updates that record instead of minting a second one that
% would report the rebuild's runtime as if it were the analysis runtime.
runRecordPath = newest_run_record_local(outputPath);
% Capture the prior analysis timings now: the saves below rewrite this file, and a
% table rebuild must not overwrite when the analysis actually ran or how long it took.
priorRecord = read_run_record_fields_local(runRecordPath, ...
    {'analysis_finished', 'runtime_seconds'});

if maxWorkers > 0
    pool = gcp('nocreate');
    if ~isempty(pool)
        assert(pool.NumWorkers == maxWorkers, ...
            'CPM:NetworkInteraction:UnexpectedExistingPool', ...
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
settings.parallel_workers_active = activeWorkers;
runRecordPath = save_analysis_run_record(outputPath, runRecordPath, runRecord);
fprintf('Parallel workers active: %d\n', activeWorkers);

checkpoint = initialize_or_load_checkpoint( ...
    checkpointPath, settings, numel(outcomeConfig));
tablesRebuiltOnly = all(~cellfun(@isempty, checkpoint.outcomes)) && ...
    ~isempty(priorRecord.analysis_finished);
if tablesRebuiltOnly
    fprintf(['All outcomes already complete; rebuilding result tables from ' ...
        'the checkpoint without rerunning any permutation.\n']);
end
runRecord.completed_result_rows = sum(~cellfun(@isempty, checkpoint.outcomes));
runRecordPath = save_analysis_run_record(outputPath, runRecordPath, runRecord);

matsFull = load_first_var_local(connectomePath);
validateattributes(matsFull, {'numeric'}, ...
    {'3d', 'real', 'finite', 'nonempty'});
assert(isequal(size(matsFull), [416, 416, 67]), ...
    'CPM:NetworkInteraction:UnexpectedConnectomeSize', ...
    'LSD_difference must be 416-by-416-by-67.');
covarsFull = load_first_var_local(covariatePath);
assert(isnumeric(covarsFull) && isequal(size(covarsFull), [67, 5]) && ...
    all(isfinite(covarsFull(:))), ...
    'CPM:NetworkInteraction:InvalidCovariates', ...
    ['covars_indicator.mat must be finite 67-by-5: ' ...
     '[study_2, study_3, sex, age, mean_FD].']);
groupingFull = load_first_var_local(groupingPath);
assert(isequal(size(groupingFull), [67, 4]), ...
    'CPM:NetworkInteraction:InvalidGrouping', ...
    'covars.mat must be 67-by-4.');
studyLabelsFull = groupingFull(:, 1);

for outcomeIndex = 1:numel(outcomeConfig)
    config = outcomeConfig(outcomeIndex);
    if ~isempty(checkpoint.outcomes{outcomeIndex})
        fprintf('\nSkipping completed outcome: %s\n', config.label);
        continue
    end

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
        'CPM:NetworkInteraction:BehaviorCountMismatch', ...
        '%s contains %d values; expected %d.', ...
        behaviorSource, numel(behavior), numel(subjectRows));

    fprintf('\n========================================\n');
    fprintf('Outcome: %s (n=%d, %s)\n', ...
        config.label, numel(subjectRows), config.correlation_type);
    fprintf('Behavior source: %s\n', behaviorSource);
    fprintf('========================================\n');

    rng(cfg.rng.cv_seed, cfg.rng.cv_generator);
    fixedPartition = cvpartition(numel(subjectRows), ...
        'KFold', NO_FOLDS);
    [behaviorColumns, freedmanLane] = ...
        freedman_lane_columns(behavior, covars, covariateNames, ...
            noIterations, cfg.rng);

    consensus = cpm_build_consensus_masks( ...
        mats, behaviorColumns, covars, config.correlation_type, ...
        fixedPartition, P_THRESHOLD, CONSENSUS_FRACTION, maxWorkers);
    assert(consensus.layout.uses_unique_edges && ...
        consensus.layout.no_edges == TOTAL_POSSIBLE_EDGES, ...
        'CPM:NetworkInteraction:UnexpectedEdgeLayout', ...
        'Expected 86,320 unique undirected edges.');

    positiveCounts = cellfun(@(x) numel(x.pos), consensus.edge_sets);
    networkPositiveCounts = count_network_pairs( ...
        consensus.edge_sets, consensus.layout, networkDefinition);
    assert(isequal(sum(double(networkPositiveCounts), 2), ...
        positiveCounts), ...
        'CPM:NetworkInteraction:NetworkCountMismatch', ...
        'Network-pair counts do not sum to positive consensus totals.');

    outcome = struct;
    outcome.behavior = config.label;
    outcome.behavior_source = behaviorSource;
    outcome.behavior_path = behaviorPath;
    outcome.behavior_sha256 = file_sha256(behaviorPath);
    outcome.behavior_fallback_used = fallbackUsed;
    outcome.no_subjects = numel(subjectRows);
    outcome.subject_rows = subjectRows;
    outcome.correlation_type = config.correlation_type;
    outcome.covariate_columns = covariateNames;
    outcome.study_reference = referenceStudy;
    outcome.observed_positive_count = positiveCounts(1);
    outcome.null_positive_counts = positiveCounts(2:end);
    outcome.network_pair_index = (1:networkDefinition.no_pairs)';
    outcome.network_pair_network1 = networkDefinition.network1;
    outcome.network_pair_network2 = networkDefinition.network2;
    outcome.network_pair_possible_edges = networkDefinition.possible_edges;
    outcome.network_positive_counts = networkPositiveCounts;
    outcome.observed_positive_mask = consensus.observed_pos_mask;
    outcome.fixed_fold_test_indices = fold_test_indices(fixedPartition);
    outcome.permutation_seeds = (201:(noIterations + 200))';
    outcome.freedman_lane = freedmanLane;

    checkpoint.outcomes{outcomeIndex} = outcome;
    checkpoint.completed_behaviors{end + 1} = config.label;
    checkpoint.last_completed_at = timestamp_now_local();
    cpm_atomic_save(checkpointPath, struct('checkpoint',checkpoint));

    fprintf('Observed positive consensus edges: %d\n', positiveCounts(1));
    fprintf('Null count: mean %.2f, median %.1f, range %d-%d\n', ...
        mean(outcome.null_positive_counts), ...
        median(outcome.null_positive_counts), ...
        min(outcome.null_positive_counts), ...
        max(outcome.null_positive_counts));

    runRecord.completed_result_rows = ...
        sum(~cellfun(@isempty, checkpoint.outcomes));
    runRecordPath = save_analysis_run_record( ...
        outputPath, runRecordPath, runRecord);
    clear consensus behaviorColumns networkPositiveCounts
end

[outcomeSummary, networkPairResults, nullCounts] = ...
    build_result_tables(checkpoint.outcomes, noIterations, ...
    TOTAL_POSSIBLE_EDGES);
writetable(outcomeSummary, summaryCsv);
writetable(networkPairResults, pairCsv);
writetable(nullCounts, nullCountCsv);
outcomeResults = checkpoint.outcomes;
analysisSettings = settings;
save(resultsMat, 'analysisSettings', 'outcomeResults', ...
    'outcomeSummary', 'networkPairResults', '-v7.3');
write_summary_text(summaryTxt, outcomeSummary, networkPairResults, settings);
build_network_interaction_workbook(outputPath);
verify_network_interaction_statistics(outputPath, noIterations);

runRecord.status = 'ANALYSIS_COMPLETE_AWAITING_POSTPROCESS';
runRecord.completed_result_rows = height(outcomeSummary);
runRecord.network_pair_result_rows = height(networkPairResults);
runRecord.null_edge_count_rows = height(nullCounts);
runRecordPath = save_analysis_run_record( ...
    outputPath, runRecordPath, runRecord);
try
    postprocessing = run_python_postprocessing( ...
        cfg.runtime.python_executable, scriptDir, outputPath, noIterations);
    runRecord.render_command = postprocessing.render_command;
    runRecord.verify_command = postprocessing.verify_command;
    runRecord.figure_verification = 'PASS';
catch exception
    runRecord.status = 'POSTPROCESSING_FAILED';
    runRecord.postprocessing_error = exception.message;
    save_analysis_run_record(outputPath, runRecordPath, runRecord);
    rethrow(exception)
end

runRecord.status = 'COMPLETE';
if tablesRebuiltOnly
    % Preserve the original analysis timings; only note the rebuild.
    runRecord.analysis_finished = priorRecord.analysis_finished;
    runRecord.runtime_seconds = priorRecord.runtime_seconds;
    runRecord.tables_rebuilt_at = timestamp_now_local();
else
    runRecord.analysis_finished = timestamp_now_local();
    runRecord.runtime_seconds = toc(runTimer);
end
runRecord.completed_result_rows = height(outcomeSummary);
runRecord.network_pair_result_rows = height(networkPairResults);
runRecord.null_edge_count_rows = height(nullCounts);
runRecord.zero_edge_null_masks_included = ...
    sum(outcomeSummary.ZeroEdgeNullMasks);
runRecordPath = save_analysis_run_record(outputPath, runRecordPath, runRecord);
save(fullfile(outputPath, 'run_record_state.mat'), ...
    'runRecord', 'runRecordPath');

fprintf('\nAnalysis, tables, and figure verification complete: %s\n', ...
    timestamp_now_local());
fprintf('Runtime: %.3f seconds\n', runRecord.runtime_seconds);
fprintf('Network-pair results: %s\n', pairCsv);
fprintf('Outcome summary: %s\n', summaryCsv);
fprintf('Null edge counts: %s\n', nullCountCsv);
fprintf('MAT results: %s\n', resultsMat);
fprintf('Publication figures: %s\n', fullfile(outputPath, 'figures'));
clear diaryCleanup
diary off
end


function details = run_python_postprocessing( ...
        pythonExecutable, scriptDir, outputPath, noIterations)
renderer = fullfile(scriptDir, 'render_network_interaction_figures.py');
verifier = fullfile(scriptDir, 'verify_network_interaction_outputs.py');
assert(isfile(renderer) && isfile(verifier), ...
    'CPM:NetworkInteraction:MissingPostprocessor', ...
    'Renderer or Python verifier is missing from the publication pipeline.');
renderCommand = sprintf('%s %s --output-dir %s', ...
    shell_quote_local(pythonExecutable), shell_quote_local(renderer), ...
    shell_quote_local(outputPath));
verifyCommand = sprintf([ ...
    '%s %s --output-dir %s --expected-permutations %d ' ...
    '--require-workbook'], ...
    shell_quote_local(pythonExecutable), shell_quote_local(verifier), ...
    shell_quote_local(outputPath), noIterations);
fprintf('Rendering publication figures: %s\n', renderCommand);
[renderStatus, renderOutput] = system(renderCommand, '-echo');
assert(renderStatus == 0, 'CPM:NetworkInteraction:RenderingFailed', ...
    'Publication figure rendering failed: %s', renderOutput);
fprintf('Verifying publication figures: %s\n', verifyCommand);
[verifyStatus, verifyOutput] = system(verifyCommand, '-echo');
assert(verifyStatus == 0, 'CPM:NetworkInteraction:FigureVerificationFailed', ...
    'Publication figure verification failed: %s', verifyOutput);
details = struct('render_command',renderCommand, ...
    'verify_command',verifyCommand);
end


function quoted = shell_quote_local(value)
value = char(string(value));
if ispc
    quoted = ['"' strrep(value, '"', '""') '"'];
else
    quoted = ['''' strrep(value, '''', '''"''"''') ''''];
end
end


function config = build_outcome_config(requestedOutcomes)
available = struct( ...
    'label', {'GDE', 'MEQ30', 'VRS'}, ...
    'preferred_behavior', {'LSD_GDE_difference', ...
        'LSD_MEQ30_difference', 'LSD_VRS_difference'}, ...
    'fallback_behavior', {'LSD_GDE', 'LSD_MEQ30', 'LSD_VRS'}, ...
    'correlation_type', {'Spearman', 'Pearson', 'Pearson'}, ...
    'use_meq30_subset', {false, true, false});
requestedOutcomes = cellstr(string(requestedOutcomes));
config = available(ismember({available.label}, requestedOutcomes));
assert(numel(config) == 3 && isequal({config.label}, requestedOutcomes) && ...
    isequal(requestedOutcomes, {'GDE','MEQ30','VRS'}), ...
    'CPM:NetworkInteraction:ConfiguredOutcomes', ...
    'Network interactions require GDE, MEQ30, and VRS in that order.');
end


function definition = load_network_definition(pathName, expectedNodes)
assert(isfile(pathName), 'CPM:NetworkInteraction:MissingAtlas', ...
    'Network atlas is missing: %s', pathName);
atlas = readtable(pathName, 'FileType', 'text', 'Delimiter', '\t', ...
    'ReadVariableNames', false);
nodeIds = atlas{:, 1};
labels = strtrim(string(atlas{:, 2}));
assert(height(atlas) == expectedNodes && ...
    isequal(nodeIds(:), (1:expectedNodes)'), ...
    'CPM:NetworkInteraction:InvalidAtlasNodes', ...
    'Atlas must contain ordered node IDs 1:%d.', expectedNodes);
names = unique(labels, 'stable');
assert(numel(names) == 8 && all(strlength(names) > 0), ...
    'CPM:NetworkInteraction:InvalidAtlasNetworks', ...
    'Expected eight nonempty network labels.');
[found, nodeNetworkId] = ismember(labels, names);
assert(all(found), 'CPM:NetworkInteraction:UnmappedAtlasNode', ...
    'Every atlas node must map to a network.');

noNetworks = numel(names);
noPairs = noNetworks * (noNetworks + 1) / 2;
pairIdMatrix = zeros(noNetworks, noNetworks, 'uint16');
network1 = strings(noPairs, 1);
network2 = strings(noPairs, 1);
possibleEdges = zeros(noPairs, 1);
nodeCounts = accumarray(nodeNetworkId, 1, [noNetworks, 1]);
pairIndex = 0;
for rowNetwork = 1:noNetworks
    for columnNetwork = 1:rowNetwork
        pairIndex = pairIndex + 1;
        pairIdMatrix(rowNetwork, columnNetwork) = pairIndex;
        pairIdMatrix(columnNetwork, rowNetwork) = pairIndex;
        network1(pairIndex) = names(columnNetwork);
        network2(pairIndex) = names(rowNetwork);
        if rowNetwork == columnNetwork
            possibleEdges(pairIndex) = ...
                nodeCounts(rowNetwork) * (nodeCounts(rowNetwork) - 1) / 2;
        else
            possibleEdges(pairIndex) = ...
                nodeCounts(rowNetwork) * nodeCounts(columnNetwork);
        end
    end
end
assert(sum(possibleEdges) == expectedNodes * (expectedNodes - 1) / 2, ...
    'CPM:NetworkInteraction:CapacityMismatch', ...
    'Network-pair capacities must sum to all undirected edges.');

definition = struct;
definition.names = names;
definition.node_network_id = nodeNetworkId;
definition.no_pairs = noPairs;
definition.pair_id_matrix = pairIdMatrix;
definition.network1 = network1;
definition.network2 = network2;
definition.possible_edges = possibleEdges;
end


function counts = count_network_pairs(edgeSets, layout, definition)
[edgeRow, edgeColumn] = ind2sub( ...
    [layout.no_node, layout.no_node], layout.upper_index);
rowNetwork = definition.node_network_id(edgeRow);
columnNetwork = definition.node_network_id(edgeColumn);
edgePairId = definition.pair_id_matrix(sub2ind( ...
    size(definition.pair_id_matrix), rowNetwork, columnNetwork));
assert(numel(edgePairId) == layout.no_edges && all(edgePairId > 0), ...
    'CPM:NetworkInteraction:EdgePairMappingFailed', ...
    'Every unique edge must map to one network pair.');

counts = zeros(numel(edgeSets), definition.no_pairs, 'uint32');
for modelIndex = 1:numel(edgeSets)
    selectedPairIds = edgePairId(double(edgeSets{modelIndex}.pos));
    if ~isempty(selectedPairIds)
        counts(modelIndex, :) = uint32(accumarray( ...
            double(selectedPairIds(:)), 1, ...
            [definition.no_pairs, 1]))';
    end
end
end


function [columns, details] = ...
        freedman_lane_columns(behavior, covars, covariateNames, ...
        noIterations, rngSettings)
behavior = behavior(:);
design = [ones(numel(behavior), 1), covars];
coefficients = lsqminnorm(design, behavior);
fitted = design * coefficients;
residuals = behavior - fitted;
columns = zeros(numel(behavior), noIterations + 1, 'like', behavior);
columns(:, 1) = behavior;
for iteration = 1:noIterations
    stream = RandStream(rngSettings.internal_permutation_generator, ...
        'Seed', iteration + rngSettings.permutation_seed_offset);
    order = randperm(stream, numel(behavior));
    columns(:, iteration + 1) = fitted + residuals(order);
end
details = struct;
details.scale = 'raw outcome';
details.design_columns = ['[intercept, ' strjoin(covariateNames, ', ') ']'];
details.design_rank = rank(design);
details.design_column_count = size(design, 2);
details.residual_mean = mean(residuals);
details.residual_sum_squares = sum(residuals .^ 2);
details.generator = rngSettings.internal_permutation_generator;
details.seed_range = [1, noIterations] + ...
    rngSettings.permutation_seed_offset;
end


function indices = fold_test_indices(partition)
indices = cell(partition.NumTestSets, 1);
for fold = 1:partition.NumTestSets
    indices{fold} = find(test(partition, fold));
end
end


function checkpoint = initialize_or_load_checkpoint(pathName, settings, noOutcomes)
if ~isfile(pathName)
    checkpoint = struct;
    checkpoint.settings = settings;
    checkpoint.outcomes = cell(1, noOutcomes);
    checkpoint.completed_behaviors = {};
    checkpoint.created_at = timestamp_now_local();
    checkpoint.last_completed_at = '';
    cpm_atomic_save(pathName, struct('checkpoint',checkpoint));
    return
end
loaded = load(pathName, 'checkpoint');
assert(isfield(loaded, 'checkpoint'), ...
    'CPM:NetworkInteraction:InvalidCheckpoint', ...
    'Checkpoint does not contain the expected variable.');
checkpoint = loaded.checkpoint;
% The guard exists to stop analysis settings or RNG conventions being mixed inside
% one checkpoint. Purely presentational fields are excluded: renaming a figure
% label must not invalidate finished permutation results.
assert(isequaln(analysis_settings_only(checkpoint.settings), ...
    analysis_settings_only(settings)), ...
    'CPM:NetworkInteraction:CheckpointSettingsMismatch', ...
    ['Existing checkpoint settings differ. Use a clean output directory; ' ...
     'never mix settings or RNG conventions.']);
checkpoint.settings = settings;
end


function settings = analysis_settings_only(settings)
%ANALYSIS_SETTINGS_ONLY Drop presentation-only fields before comparing settings.
displayOnlyFields = {'display_effect_size', ...
    'parallel_workers_requested', 'parallel_workers_active'};
for index = 1:numel(displayOnlyFields)
    if isfield(settings, displayOnlyFields{index})
        settings = rmfield(settings, displayOnlyFields{index});
    end
end
end


function [summaryTable, pairTable, nullTable] = ...
        build_result_tables(outcomes, noIterations, totalPossibleEdges)
noOutcomes = numel(outcomes);
summaryParts = cell(noOutcomes, 1);
pairParts = cell(noOutcomes, 1);
nullParts = cell(noOutcomes, 1);

for outcomeIndex = 1:noOutcomes
    current = outcomes{outcomeIndex};
    assert(~isempty(current), ...
        'CPM:NetworkInteraction:IncompleteCheckpoint', ...
        'All outcomes must finish before final tables are built.');
    counts = double(current.network_positive_counts);
    totals = sum(counts, 2);
    possible = current.network_pair_possible_edges(:)';
    pairDensity = counts ./ possible;
    overallDensity = totals / totalPossibleEdges;
    enrichment = pairDensity ./ overallDensity;
    observedEnrichment = enrichment(1, :);
    nullCounts = counts(2:end, :);
    countP = (1 + sum(nullCounts >= counts(1, :), 1))' ...
        / (noIterations + 1);
    countQ = benjamini_hochberg(countP);

    behavior = repmat(string(current.behavior), numel(possible), 1);
    pairIndex = current.network_pair_index;
    network1 = paper_network_names(current.network_pair_network1);
    network2 = paper_network_names(current.network_pair_network2);
    networkPair = network1 + "-" + network2;
    observedCount = counts(1, :)';
    possibleEdges = possible';
    networkPairDensity = pairDensity(1, :)';
    observedExpectedRatio = observedEnrichment';
    log2ObservedExpectedRatio = log2(observedExpectedRatio);
    countPermutationP = countP;
    countFdrQ = countQ;
    countFdrSignificant = countQ < 0.05;
    sparseObservedCount = observedCount < 5;
    pairParts{outcomeIndex} = table(behavior, pairIndex, network1, ...
        network2, networkPair, observedCount, possibleEdges, ...
        networkPairDensity, observedExpectedRatio, ...
        log2ObservedExpectedRatio, ...
        countPermutationP, countFdrQ, countFdrSignificant, ...
        sparseObservedCount, ...
        'VariableNames', {'Outcome', 'PairIndex', 'Network1', 'Network2', ...
        'NetworkPair', 'ObservedConsensusEdges', 'PossibleEdges', ...
        'NetworkPairDensity', 'ObservedExpectedRatio', ...
        'Log2ObservedExpectedRatio', ...
        'CountPermutationP', 'CountFDRQ', 'CountFDRSignificant', ...
        'SparseObservedCount'});

    nullValues = totals(2:end);
    totalMaskPermutationP = ...
        (1 + sum(nullValues >= totals(1))) / (noIterations + 1);
    zeroEdgeNullMasks = sum(nullValues == 0);
    summaryParts{outcomeIndex} = table(string(current.behavior), ...
        current.no_subjects, string(current.correlation_type), totals(1), ...
        mean(nullValues), median(nullValues), ...
        prctile(nullValues, 2.5), prctile(nullValues, 97.5), ...
        min(nullValues), max(nullValues), totalMaskPermutationP, ...
        zeroEdgeNullMasks, sum(countFdrSignificant), ...
        'VariableNames', {'Outcome', 'NSubjects', 'CorrelationType', ...
        'ObservedPositiveConsensusEdges', 'NullMean', 'NullMedian', ...
        'NullP025', 'NullP975', 'NullMinimum', 'NullMaximum', ...
        'TotalMaskPermutationP', 'ZeroEdgeNullMasks', ...
        'CountFDRSignificantNetworkPairs'});

    nullParts{outcomeIndex} = table( ...
        repmat(string(current.behavior), noIterations, 1), ...
        (1:noIterations)', totals(2:end), ...
        'VariableNames', {'Outcome', 'Permutation', ...
        'PositiveConsensusEdgeCount'});
end
summaryTable = vertcat(summaryParts{:});
summaryTable.TotalMaskFDRQ = benjamini_hochberg( ...
    summaryTable.TotalMaskPermutationP);
summaryTable.TotalMaskFDRSignificant = ...
    summaryTable.TotalMaskFDRQ < 0.05;
pairTable = vertcat(pairParts{:});
nullTable = vertcat(nullParts{:});
end


function names = paper_network_names(names)
names = string(names);
names(names == "ASM") = "SMN";
names(names == "Subcortical") = "SBC";
end


function write_summary_text(pathName, summary, pairs, settings)
fileId = fopen(pathName, 'wt');
assert(fileId ~= -1, 'CPM:NetworkInteraction:CannotWriteSummary', ...
    'Could not create summary file: %s', pathName);
cleanup = onCleanup(@() fclose(fileId));
fprintf(fileId, 'CPM network-pair enrichment with Freedman-Lane null\n');
fprintf(fileId, 'Completed: %s\n\n', timestamp_now_local());
fprintf(fileId, 'Positive LSD-PCB difference consensus masks only.\n');
fprintf(fileId, ['Each outcome used %d behavioral residual permutations, ' ...
    '10 fixed folds, p < %.3g, and >=8/10-fold consensus.\n'], ...
    settings.no_iterations, settings.edge_p_threshold);
fprintf(fileId, ['Covariates: study 2 indicator, study 3 indicator, sex, ' ...
    'age, and mean framewise displacement.\n']);
fprintf(fileId, ['The observed / expected ratio is the displayed effect ' ...
    'size, where 1 is the share an even spread across all connections would ' ...
    'give a pair. Primary ' ...
    'network-pair inference uses the raw consensus-edge count.\n']);
fprintf(fileId, ['Raw-count pair tests are one-sided upper-tail tests; all ' ...
    'null masks are included, and Benjamini-Hochberg FDR is applied ' ...
    'across 36 pairs per outcome.\n']);
fprintf(fileId, ['Total-mask edge-count p-values are separate one-sided ' ...
    'permutation tests; their three-outcome family is BH-corrected ' ...
    'separately from the pairwise tests.\n\n']);
for row = 1:height(summary)
    fprintf(fileId, ['%s (n=%d, %s): observed %d edges; null median %.1f ' ...
        '[2.5th-97.5th percentile %.1f-%.1f]; total-mask p=%.6f, ' ...
        'q=%.6f; ' ...
        'count-FDR-significant pairs=%d.\n'], ...
        summary.Outcome(row), summary.NSubjects(row), ...
        summary.CorrelationType(row), ...
        summary.ObservedPositiveConsensusEdges(row), ...
        summary.NullMedian(row), summary.NullP025(row), ...
        summary.NullP975(row), summary.TotalMaskPermutationP(row), ...
        summary.TotalMaskFDRQ(row), ...
        summary.CountFDRSignificantNetworkPairs(row));
    subset = pairs(pairs.Outcome == summary.Outcome(row) & ...
        pairs.CountFDRSignificant, :);
    if isempty(subset)
        fprintf(fileId, ...
            '  No network pair survived raw-count within-outcome FDR.\n');
    else
        for pairRow = 1:height(subset)
            fprintf(fileId, ...
                '  %s: %.3fx, count p=%.6f, count q=%.6f, edges=%d.\n', ...
                subset.NetworkPair(pairRow), ...
                subset.ObservedExpectedRatio(pairRow), ...
                subset.CountPermutationP(pairRow), ...
                subset.CountFDRQ(pairRow), ...
                subset.ObservedConsensusEdges(pairRow));
        end
    end
end
clear cleanup
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
        'CPM:NetworkInteraction:MissingBehavior', ...
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
assert(isfile(pathName), 'CPM:NetworkInteraction:MissingInput', ...
    'Required input is missing: %s', pathName);
loaded = load(pathName);
names = fieldnames(loaded);
assert(isscalar(names), 'CPM:NetworkInteraction:AmbiguousInput', ...
    'Expected exactly one variable in %s.', pathName);
value = loaded.(names{1});
end


function recordPath = newest_run_record_local(outputPath)
%NEWEST_RUN_RECORD_LOCAL Path of the existing run record, or '' if none.
listing = dir(fullfile(outputPath, 'run_record_*.txt'));
if isempty(listing)
    recordPath = '';
    return
end
[~, newest] = max([listing.datenum]);
recordPath = fullfile(outputPath, listing(newest).name);
end


function fields = read_run_record_fields_local(recordPath, names)
%READ_RUN_RECORD_FIELDS_LOCAL Read named key=value pairs from a run record.
fields = struct('analysis_finished', '', 'runtime_seconds', 0);
if isempty(recordPath) || ~isfile(recordPath)
    return
end
text = fileread(recordPath);
lines = strsplit(text, newline);
for index = 1:numel(names)
    key = names{index};
    match = lines(startsWith(lines, [key '=']));
    if isempty(match)
        continue
    end
    value = extractAfter(match{end}, [key '=']);
    numeric = str2double(value);
    if ~isnan(numeric)
        fields.(key) = numeric;
    else
        fields.(key) = char(value);
    end
end
end


function value = timestamp_now_local()
value = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS'));
end
