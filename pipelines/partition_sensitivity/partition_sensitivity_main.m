function results = partition_sensitivity_main(cfg, options)
%PARTITION_SENSITIVITY_MAIN Repeat the primary model across 100 partitions.
%
% This analysis produces the data for Supplementary Figure S1. It uses
% the positive covariate-adjusted LSD-placebo difference model for all six
% outcomes and CV seeds 123:222 by default.

arguments
    cfg (1,1) struct
    options.no_iterations (1,1) double = NaN
    options.workers (1,1) double = NaN
    options.output_path (1,1) string = ""
    options.seeds (1,:) double = []
end

pipelineDir = fileparts(mfilename('fullpath'));
repositoryRoot = fileparts(fileparts(pipelineDir));
addpath(pipelineDir, fullfile(repositoryRoot, 'scripts', 'common'), ...
    fullfile(repositoryRoot, 'scripts', 'utilities'));
if isnan(options.no_iterations), noIterations = cfg.partition.iterations;
else, noIterations = options.no_iterations; end
if isnan(options.workers), workers = cfg.runtime.workers;
else, workers = options.workers; end
if isempty(options.seeds), seeds = cfg.partition.seeds;
else, seeds = options.seeds; end
if options.output_path == ""
    outputPath = fullfile(cfg.output_root, 'partition_sensitivity');
else
    outputPath = char(options.output_path);
end
if ~isfolder(outputPath), mkdir(outputPath); end
checkpointPath = fullfile(outputPath, 'partition_sensitivity_checkpoint.mat');

[inputs, code] = cpm_analysis_dependencies(cfg, {pipelineDir}, ...
    {fullfile(cfg.input_root,'connectomes','LSD_difference.mat')});
settings = struct('outcomes',cfg.outcomes, 'selected', {cfg.partition.outcomes}, ...
    'seeds',seeds, 'iterations',noIterations, 'rng',cfg.rng, ...
    'folds',cfg.model.folds,'threshold',cfg.model.edge_threshold);
cpm_checkpoint_guard(outputPath, settings, inputs, code);
started = tic;
validateattributes(noIterations, {'numeric'}, {'scalar','integer','>=',2});
validateattributes(workers, {'numeric'}, {'scalar','integer','>=',0});
validateattributes(seeds, {'numeric'}, {'vector','integer','nonnegative'});
if workers > 0
    pool = gcp('nocreate');
    if ~isempty(pool) && pool.NumWorkers ~= workers, delete(pool); pool = []; end
    if isempty(pool), parpool('Processes', workers); end
end

record = struct('status','RUNNING', 'analysis','partition_sensitivity', ...
    'entry_script',[mfilename('fullpath') '.m'], 'input_path',cfg.input_root, ...
    'output_path',outputPath, 'cv_seeds',mat2str(seeds), ...
    'no_iterations',noIterations, 'parallel_workers_requested',workers, ...
    'randomization_scheme', sprintf([ ...
        'CV rng(seed,%s); permutations rng(iteration+%d,%s)'], ...
        cfg.rng.cv_generator, cfg.rng.permutation_seed_offset, ...
        cfg.rng.internal_permutation_generator), ...
    'checkpoint_mat',checkpointPath, 'completed_result_rows',0);
record.invocation_started = char(datetime('now'));
record.covariates_input_file = cfg.model.covariate_file;
record.covariate_columns = cfg.model.covariate_columns;
record.parallel_workers_active = workers;
record.behaviors = cfg.partition.outcomes;
record.folds = cfg.model.folds;
record.threshold = cfg.model.edge_threshold;
record.subject_subsets = 'MEQ30 20:67; other outcomes 1:67';
record.difference_behavior_fallback = 'prefer difference; standard is equivalent when placebo is zero';
record.canonical_result_files = {checkpointPath, ...
    fullfile(outputPath,'partition_sensitivity_results.csv'), ...
    fullfile(outputPath,'partition_sensitivity_results.mat')};
runRecordPath = save_analysis_run_record(outputPath, '', record);

matsFull = load_single_numeric(fullfile(cfg.input_root, 'connectomes', ...
    'LSD_difference.mat'));
covarsFull = load_single_numeric(cfg.model.covariate_file);
grouping = load_single_numeric(cfg.model.grouping_file);
studyLabelsFull = grouping(:, 1);
outcomes = cfg.outcomes(ismember({cfg.outcomes.label}, cfg.partition.outcomes));
if isfile(checkpointPath)
    loadedCheckpoint = load(checkpointPath, 'rows', 'checkpointSeeds');
    assert(isequal(loadedCheckpoint.checkpointSeeds, seeds), ...
        'CPM:Partition:CheckpointSettings', ...
        'Existing checkpoint uses different CV seeds.');
    rows = loadedCheckpoint.rows;
else
    rows = cell(0, 24);
end
for outcomeIndex = 1:numel(outcomes)
    outcome = outcomes(outcomeIndex);
    subjectRows = outcome.subject_rows;
    mats = matsFull(:, :, subjectRows);
    [covars, covariateNames, referenceStudy] = ...
        cpm_select_covariates(covarsFull(subjectRows, :), ...
            studyLabelsFull(subjectRows));
    behaviorStem = outcome.difference_behavior;
    behaviorPath = fullfile(cfg.input_root, 'behav', [behaviorStem '.mat']);
    if ~isfile(behaviorPath)
        behaviorStem = outcome.standard_behavior;
        behaviorPath = fullfile(cfg.input_root, 'behav', [behaviorStem '.mat']);
    end
    behavior = load_single_numeric(behaviorPath);
    behavior = behavior(:);
    if numel(behavior) == 67 && numel(subjectRows) < 67
        behavior = behavior(subjectRows);
    end

    for seed = seeds
        if ~isempty(rows) && any(strcmp(rows(:,1), outcome.label) & ...
                cell2mat(rows(:,2)) == seed)
            continue
        end
        [engineRows, ~] = permutation_test_cv_partition( ...
            'LSD_difference', behaviorStem, mats, behavior, ...
            cfg.model.folds, cfg.model.edge_threshold, covars, '', {}, ...
            noIterations, seed, false, outcome.correlation_type, workers, ...
            cfg.rng);
        row = engineRows(1, :);
        rows(end + 1, :) = [ ... %#ok<AGROW>
            {outcome.label, seed, numel(subjectRows), ...
             outcome.correlation_type, strjoin(covariateNames, ', '), ...
             referenceStudy}, row, ...
            {'LSD-PCB difference','With covariates','Positive'}];
        checkpointSeeds = seeds; %#ok<NASGU>
        cpm_atomic_save(checkpointPath, struct('rows',{rows}, ...
            'checkpointSeeds',checkpointSeeds));
        record.completed_result_rows = size(rows,1);
        save_analysis_run_record(outputPath, runRecordPath, record);
    end
end

names = {'Outcome','CVSeed','N','CorrelationType','CovariateColumns', ...
    'StudyReference','Analysis','K', ...
    'Threshold','CombinedR','CombinedP','PositiveR','PositiveP', ...
    'NegativeR','NegativeP','CombinedMSE','CombinedMSEP','PositiveMSE', ...
    'PositiveMSEP','NegativeMSE','NegativeMSEP','ModelForm', ...
    'Covariates','Network'};
% The engine contributes 15 columns, giving 24 columns in total.
results = cell2table(rows, 'VariableNames', names);
writetable(results, fullfile(outputPath, 'partition_sensitivity_results.csv'));
save(fullfile(outputPath, 'partition_sensitivity_results.mat'), ...
    'results', 'cfg', 'seeds', 'noIterations', '-v7.3');
record.status = 'COMPLETE';
record.runtime_seconds = toc(started);
record.invocation_finished = char(datetime('now'));
record.covariate_columns = 'per-outcome CovariateColumns and StudyReference in canonical CSV';
record.completed_result_rows = height(results);
save_analysis_run_record(outputPath, runRecordPath, record);
end

function value = load_single_numeric(pathName)
assert(isfile(pathName), 'CPM:Partition:MissingInput', ...
    'Missing input: %s', pathName);
loaded = load(pathName);
names = fieldnames(loaded);
assert(isscalar(names) && isnumeric(loaded.(names{1})), ...
    'CPM:Partition:InputFormat', ...
    '%s must contain one numeric variable.', pathName);
value = loaded.(names{1});
assert(all(isfinite(value), 'all'), 'CPM:Partition:NonfiniteInput', ...
    '%s contains NaN or Inf.', pathName);
end
