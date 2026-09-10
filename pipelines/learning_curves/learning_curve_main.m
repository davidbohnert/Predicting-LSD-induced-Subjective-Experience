function summary = learning_curve_main(cfg, options)
%LEARNING_CURVE_MAIN Learning curve for the paper's primary model.
%
% Sample-size dependence, restricted to the designated
% primary analysis: the LSD-PCB difference connectome, adjusted for
% indicator-coded covariates. All six experience dimensions are included,
% in the Table 1 order. The engine retains combined, positive, and negative
% networks; the positive network is the reported one.

arguments
    cfg (1,1) struct
    options.no_iterations (1,1) double = NaN
    options.workers (1,1) double = NaN
    options.output_path (1,1) string = ""
    options.repetitions (1,1) double = NaN
    options.standard_sample_sizes (1,:) double = []
    options.meq30_sample_sizes (1,:) double = []
end
scriptRoot = fileparts(mfilename('fullpath'));
repositoryRoot = fileparts(fileparts(scriptRoot));
if options.output_path == ""
    runRoot = fullfile(cfg.output_root, 'learning_curves');
else
    runRoot = char(options.output_path);
end
outputPath = fullfile(runRoot, 'output');
inputPath = cfg.input_root;
engineFile = fullfile(scriptRoot, 'permutation_test_cv_sample_size.m');
utilitiesPath = fullfile(repositoryRoot, 'scripts', 'utilities');
commonPath = fullfile(repositoryRoot, 'scripts', 'common');
addpath(scriptRoot, utilitiesPath, commonPath);
if ~isfolder(outputPath), mkdir(outputPath); end

rawCsv = fullfile(outputPath, 'primary_learning_curve_raw_results.csv');
summaryCsv = fullfile(outputPath, 'primary_learning_curve_summary.csv');
checkpointMat = fullfile(outputPath, ...
    'primary_learning_curve_checkpoint.mat');
runLog = fullfile(outputPath, 'primary_learning_curve_run_log.txt');
planPath = fullfile(runRoot, 'primary_learning_curve_plan.mat');

K_FOLDS = cfg.model.folds;
EDGE_THRESHOLD = cfg.model.edge_threshold;
if isnan(options.no_iterations), NO_ITERATIONS = cfg.learning.iterations;
else, NO_ITERATIONS = options.no_iterations; end
if isnan(options.workers), POOL_WORKERS = cfg.runtime.workers;
else, POOL_WORKERS = options.workers; end
if isnan(options.repetitions), NO_REPETITIONS = cfg.learning.repetitions;
else, NO_REPETITIONS = options.repetitions; end
CHECKPOINT_EVERY = 5;
MEQ30_SUBSET = 20:67;
% Table 1 order (rows ordered by effect size in the manuscript).
OUTCOMES = {'GDE', 'MEQ30', 'VRS', 'OBN', 'BDE', 'AED'};
if isempty(options.standard_sample_sizes)
    STANDARD_SAMPLE_SIZES = cfg.learning.standard_sample_sizes;
else
    STANDARD_SAMPLE_SIZES = options.standard_sample_sizes;
end
if isempty(options.meq30_sample_sizes)
    MEQ30_SAMPLE_SIZES = cfg.learning.meq30_sample_sizes;
else
    MEQ30_SAMPLE_SIZES = options.meq30_sample_sizes;
end
CONNECTOME_FORM = 'LSD_difference';
COVARIATE_SETTING = 'With covariates';
EXPECTED_ROWS = 5 * numel(STANDARD_SAMPLE_SIZES) * NO_REPETITIONS + ...
    numel(MEQ30_SAMPLE_SIZES) * NO_REPETITIONS;
validateattributes(NO_REPETITIONS, {'numeric'}, ...
    {'scalar','integer','positive'});

assert(isfile(engineFile), 'CPM:PrimaryCurve:EngineMissing', ...
    'Missing learning-curve engine: %s', engineFile);

[inputs, code] = cpm_analysis_dependencies(cfg, {scriptRoot}, ...
    {fullfile(inputPath,'connectomes','LSD_difference.mat')});
settings = struct('outcomes',{OUTCOMES}, 'iterations',NO_ITERATIONS, ...
    'repetitions',NO_REPETITIONS, 'sizes',STANDARD_SAMPLE_SIZES, ...
    'meq_sizes',MEQ30_SAMPLE_SIZES, 'rng',cfg.rng, ...
    'folds',K_FOLDS, 'threshold',EDGE_THRESHOLD, 'meq_rows',MEQ30_SUBSET);
cpm_checkpoint_guard(outputPath, settings, inputs, code);
diary(runLog);
diary on;
diaryCleanup = onCleanup(@() diary('off'));
fprintf('Primary-model learning-curve analysis started: %s\n', ...
    timestamp_now());
fprintf(['Scope: six outcomes; %s connectome; %s only; k=%d; ' ...
    'threshold=%.3f; iterations=%d; repetitions=%d; workers=%d\n'], ...
    CONNECTOME_FORM, lower(COVARIATE_SETTING), K_FOLDS, ...
    EDGE_THRESHOLD, NO_ITERATIONS, NO_REPETITIONS, POOL_WORKERS);
fprintf('Outcome order (Table 1): %s\n', strjoin(OUTCOMES, ', '));
fprintf('Standard sample sizes: %s\n', mat2str(STANDARD_SAMPLE_SIZES));
fprintf('MEQ30 sample sizes: %s (established maximum N=48)\n', ...
    mat2str(MEQ30_SAMPLE_SIZES));

% Indicator-coded covariates are mandatory for every adjusted analysis.
[covarsIndicator, covarsVariable] = load_one_numeric( ...
    cfg.model.covariate_file);
assert(isequal(size(covarsIndicator), [67 5]) && ...
    all(isfinite(covarsIndicator(:))), 'CPM:PrimaryCurve:Covariates', ...
    'covars_indicator.mat must contain a finite 67-by-5 matrix.');
assert(all(ismember(covarsIndicator(:, 1:3), [0 1]), 'all'), ...
    'CPM:PrimaryCurve:CovariateCoding', ...
    'Columns 1-3 of covars_indicator.mat must be binary.');

% Study labels are used only to build stratified folds, never as a regressor.
[covarsGrouping, ~] = load_one_numeric(fullfile(inputPath, 'covars.mat'));
assert(isequal(size(covarsGrouping), [67 4]), ...
    'CPM:PrimaryCurve:GroupingFile', ...
    'covars.mat must contain a 67-by-4 matrix.');
studyLabels = covarsGrouping(:, 1);

plan = learning_curve_plan(planPath, studyLabels, ...
    studyLabels(MEQ30_SUBSET), STANDARD_SAMPLE_SIZES, ...
    MEQ30_SAMPLE_SIZES, NO_REPETITIONS, cfg.rng);

[differenceMatsFull, differenceVariable] = load_one_numeric(fullfile( ...
    inputPath, 'connectomes', 'LSD_difference.mat'));
validate_connectome(differenceMatsFull, 'LSD_difference', 67);
fprintf('Loaded connectome: LSD_difference=%s\n', differenceVariable);
fprintf('Loaded covariates: %s (%s), columns = ', ...
    covarsVariable, mat2str(size(covarsIndicator)));
fprintf('[study_2, study_3, sex, age, mean_FD]\n');

entryScript = [mfilename('fullpath') '.m'];
runRecord = make_run_record(entryScript, engineFile, inputPath, ...
    outputPath, rawCsv, summaryCsv, checkpointMat, runLog, planPath, ...
    OUTCOMES, STANDARD_SAMPLE_SIZES, MEQ30_SAMPLE_SIZES, ...
    MEQ30_SUBSET, K_FOLDS, EDGE_THRESHOLD, NO_ITERATIONS, POOL_WORKERS, ...
    NO_REPETITIONS, EXPECTED_ROWS, CONNECTOME_FORM, COVARIATE_SETTING, ...
    plan, cfg.rng);
runRecordPath = save_analysis_run_record(outputPath, '', runRecord);

rawResults = load_checkpoint(checkpointMat, empty_results_table());
assert(height(rawResults) <= EXPECTED_ROWS, ...
    'CPM:PrimaryCurve:CheckpointTooLarge', ...
    'Checkpoint contains more than the expected %d rows.', EXPECTED_ROWS);
runRecord.completed_result_rows = height(rawResults);
runRecordPath = save_analysis_run_record( ...
    outputPath, runRecordPath, runRecord);

activeWorkers = ensure_pool(POOL_WORKERS);
runRecord.parallel_workers_active = activeWorkers;
runRecordPath = save_analysis_run_record( ...
    outputPath, runRecordPath, runRecord);

analysisStart = tic;
for outcomeIndex = 1:numel(OUTCOMES)
    outcome = OUTCOMES{outcomeIndex};
    [cohortPlan, differenceMats, cohortCovars, differenceBehavior, ...
        differenceBehaviorName, cohortLabel] = prepare_outcome( ...
        outcome, inputPath, plan, differenceMatsFull, ...
        covarsIndicator, studyLabels, MEQ30_SUBSET);

    fprintf('\n============================================================\n');
    fprintf('Outcome %d/%d: %s (%s)\n', outcomeIndex, ...
        numel(OUTCOMES), outcome, cohortLabel);
    fprintf('Behavior vector: %s\n', differenceBehaviorName);
    fprintf('Sample sizes: %s\n', mat2str(cohortPlan.sample_sizes));

    for sizeIndex = 1:numel(cohortPlan.sample_sizes)
        sampleN = cohortPlan.sample_sizes(sizeIndex);
        studyCounts = cohortPlan.retained_study_counts(sizeIndex, :);
        fprintf('\n---- %s N=%d; studies=%s ----\n', ...
            outcome, sampleN, mat2str(studyCounts));

        for repetition = 1:NO_REPETITIONS
            if result_exists(rawResults, outcome, sampleN, repetition)
                continue
            end
            selected = cohortPlan.sample_indices{sizeIndex, repetition};
            sampleSeed = cohortPlan.sample_seeds(repetition);
            cvSeed = plan.cv_seeds(repetition);

            mats = differenceMats(:, :, selected);
            behavior = differenceBehavior(selected);
            modelCovars = cohortCovars(selected, :);

            fprintf('%s N=%d rep=%d/%d sampleSeed=%d cvSeed=%d\n', ...
                outcome, sampleN, repetition, NO_REPETITIONS, ...
                sampleSeed, cvSeed);
            engineRows = permutation_test_cv_sample_size( ...
                CONNECTOME_FORM, differenceBehaviorName, mats, behavior, ...
                K_FOLDS, EDGE_THRESHOLD, modelCovars, '', {}, ...
                NO_ITERATIONS, cvSeed, false, POOL_WORKERS, cfg.rng);
            assert(isequal(size(engineRows), [1 15]), ...
                'CPM:PrimaryCurve:EngineShape', ...
                'The engine returned an unexpected row.');
            rawResults(end+1, :) = result_row(outcome, cohortLabel, ...
                CONNECTOME_FORM, differenceBehaviorName, ...
                COVARIATE_SETTING, sampleN, repetition, sampleSeed, ...
                cvSeed, studyCounts, K_FOLDS, EDGE_THRESHOLD, ...
                engineRows); %#ok<AGROW>

            if mod(repetition, CHECKPOINT_EVERY) == 0 || ...
                    repetition == NO_REPETITIONS
                checkpoint_results(rawResults, rawCsv, checkpointMat, plan);
                runRecord.completed_result_rows = height(rawResults);
                runRecordPath = save_analysis_run_record( ...
                    outputPath, runRecordPath, runRecord);
                fprintf('CHECKPOINT: %d/%d rows complete (%.1f min elapsed).\n', ...
                    height(rawResults), EXPECTED_ROWS, toc(analysisStart)/60);
            end
        end
    end
end

assert(height(rawResults) == EXPECTED_ROWS, ...
    'CPM:PrimaryCurve:IncompleteResults', ...
    'Expected %d rows but found %d.', EXPECTED_ROWS, height(rawResults));
assert(height(unique(rawResults(:, ...
    {'Outcome', 'SampleN', 'Repetition'}))) == EXPECTED_ROWS, ...
    'CPM:PrimaryCurve:DuplicateResults', ...
    'The completed checkpoint contains duplicate result keys.');

summary = summarize_learning_curves(runRoot);
render_learning_curve_figure(runRoot);
runRecord.status = 'COMPLETE';
runRecord.invocation_finished = timestamp_now();
runRecord.runtime_seconds = toc(analysisStart);
runRecord.completed_result_rows = height(rawResults);
runRecord.completed_summary_rows = height(summary);
runRecord.results_text_sha256 = file_sha256(rawCsv);
runRecord.results_summary_sha256 = file_sha256(summaryCsv);
save_analysis_run_record(outputPath, runRecordPath, runRecord);

fprintf('\nPrimary-model learning-curve analysis COMPLETE.\n');
fprintf('Raw result rows: %d\n', height(rawResults));
fprintf('Summary rows: %d\n', height(summary));
fprintf('Analysis runtime: %.1f seconds (%.2f hours)\n', ...
    runRecord.runtime_seconds, runRecord.runtime_seconds / 3600);
fprintf('Raw results: %s\n', rawCsv);
fprintf('Summary: %s\n', summaryCsv);
fprintf('Run record: %s\n', runRecordPath);
clear diaryCleanup
diary('off');
end


function [cohortPlan, differenceMats, cohortCovars, differenceBehavior, ...
        differenceBehaviorName, cohortLabel] = prepare_outcome( ...
        outcome, inputPath, plan, differenceMatsFull, covarsIndicator, ...
        studyLabels, meq30Subset)
standardBehaviorName = ['LSD_' outcome];
[differenceBehavior, differenceBehaviorName] = ...
    load_behavior_with_fallback(inputPath, ...
    [standardBehaviorName '_difference'], standardBehaviorName);

if strcmp(outcome, 'MEQ30')
    cohortPlan = plan.meq30;
    differenceMats = differenceMatsFull(:, :, meq30Subset);
    cohortCovars = cpm_select_covariates( ...
        covarsIndicator(meq30Subset, :), studyLabels(meq30Subset));
    % MEQ30 behavior files are already stored as the 48-row subset.
    expectedN = 48;
    cohortLabel = 'MEQ30 subset (LSD rows 20:67)';
else
    cohortPlan = plan.standard;
    differenceMats = differenceMatsFull;
    cohortCovars = cpm_select_covariates(covarsIndicator, studyLabels);
    expectedN = 67;
    cohortLabel = 'Full LSD cohort';
end

assert(numel(differenceBehavior) == expectedN && ...
    all(isfinite(differenceBehavior)), ...
    'CPM:PrimaryCurve:BehaviorLength', ...
    '%s behavior input must contain %d finite values.', ...
    outcome, expectedN);
end


function [behavior, sourceName] = load_behavior_with_fallback( ...
        inputPath, preferredName, fallbackName)
% A missing *_difference file means every placebo rating on that scale was
% zero, so the standard vector already is the difference vector.
preferredPath = fullfile(inputPath, 'behav', [preferredName '.mat']);
if isfile(preferredPath)
    sourceName = preferredName;
else
    sourceName = fallbackName;
end
[behavior, ~] = load_one_numeric(fullfile( ...
    inputPath, 'behav', [sourceName '.mat']));
behavior = behavior(:);
end


function row = result_row(outcome, cohort, connectomeForm, behaviorForm, ...
        covariateSetting, sampleN, repetition, sampleSeed, cvSeed, ...
        studyCounts, kFolds, edgeThreshold, engineRows)
row = {string(outcome), string(cohort), string(connectomeForm), ...
    string(behaviorForm), string(covariateSetting), sampleN, repetition, ...
    sampleSeed, cvSeed, studyCounts(1), studyCounts(2), studyCounts(3), ...
    kFolds, edgeThreshold, engineRows{4}, engineRows{5}, ...
    engineRows{6}, engineRows{7}, engineRows{8}, engineRows{9}, ...
    engineRows{10}, engineRows{11}, engineRows{12}, engineRows{13}, ...
    engineRows{14}, engineRows{15}};
end


function tf = result_exists(results, outcome, sampleN, repetition)
tf = any(results.Outcome == string(outcome) & ...
    results.SampleN == sampleN & results.Repetition == repetition);
end


function tableValue = empty_results_table
variableTypes = [repmat({'string'}, 1, 5), repmat({'double'}, 1, 21)];
variableNames = {'Outcome','Cohort','ConnectomeForm','BehaviorFormUsed', ...
    'CovariateSetting','SampleN','Repetition','SampleSeed','CVSeed', ...
    'Study1N','Study2N','Study3N','KFolds','EdgeThreshold', ...
    'RCombined','PCombined','RPositive','PPositive','RNegative','PNegative', ...
    'MSECombined','PMSECombined','MSEPositive','PMSEPositive', ...
    'MSENegative','PMSENegative'};
tableValue = table('Size', [0 numel(variableNames)], ...
    'VariableTypes', variableTypes, 'VariableNames', variableNames);
end


function results = load_checkpoint(checkpointMat, emptyTable)
if isfile(checkpointMat)
    loaded = load(checkpointMat, 'rawResults');
    results = loaded.rawResults;
    assert(isequal(results.Properties.VariableNames, ...
        emptyTable.Properties.VariableNames), ...
        'CPM:PrimaryCurve:CheckpointHeader', ...
        'The existing checkpoint has unexpected columns.');
else
    results = emptyTable;
end
end


function checkpoint_results(rawResults, rawCsv, checkpointMat, plan)
temporaryCsv = [tempname(fileparts(rawCsv)) '.csv'];
writetable(rawResults, temporaryCsv);
[ok, message] = movefile(temporaryCsv, rawCsv, 'f');
assert(ok, 'CPM:PrimaryCurve:CsvCheckpoint', ...
    'Could not checkpoint %s: %s', rawCsv, message);
temporaryMat = [tempname(fileparts(checkpointMat)) '.mat'];
save(temporaryMat, 'rawResults', 'plan', '-v7.3');
[ok, message] = movefile(temporaryMat, checkpointMat, 'f');
assert(ok, 'CPM:PrimaryCurve:MatCheckpoint', ...
    'Could not checkpoint %s: %s', checkpointMat, message);
end


function activeWorkers = ensure_pool(workerCount)
pool = gcp('nocreate');
if workerCount == 0
    activeWorkers = 0;
    return
end
if ~isempty(pool) && pool.NumWorkers ~= workerCount
    delete(pool);
    pool = [];
end
if isempty(pool), pool = parpool('Processes', workerCount); end
assert(pool.NumWorkers == workerCount, ...
    'CPM:PrimaryCurve:WorkerCount', ...
    'Expected %d workers; started %d.', workerCount, pool.NumWorkers);
activeWorkers = pool.NumWorkers;
end


function validate_connectome(value, name, expectedSubjects)
assert(isnumeric(value) && ndims(value) == 3 && ...
    size(value, 1) == 416 && size(value, 2) == 416 && ...
    size(value, 3) == expectedSubjects && all(isfinite(value(:))), ...
    'CPM:PrimaryCurve:InvalidConnectome', ...
    '%s must be 416-by-416-by-%d and finite.', name, expectedSubjects);
end


function [value, variableName] = load_one_numeric(pathName)
assert(isfile(pathName), 'CPM:PrimaryCurve:MissingInput', ...
    'Missing input: %s', pathName);
loaded = load(pathName);
names = fieldnames(loaded);
assert(isscalar(names), 'CPM:PrimaryCurve:AmbiguousInput', ...
    '%s must contain exactly one variable.', pathName);
variableName = names{1};
value = loaded.(variableName);
assert(isnumeric(value), 'CPM:PrimaryCurve:NonNumericInput', ...
    '%s must contain numeric data.', pathName);
end


function record = make_run_record(entryScript, engineFile, inputPath, ...
        outputPath, rawCsv, summaryCsv, checkpointMat, runLog, planPath, ...
        outcomes, standardSizes, meqSizes, meqSubset, ...
        kFolds, edgeThreshold, noIterations, workers, repetitions, ...
        expectedRows, connectomeForm, covariateSetting, plan, rngSettings)
record = struct;
record.status = 'RUNNING';
record.invocation_started = timestamp_now();
record.analysis = 'primary_model_internal_validation_learning_curve';
record.purpose = [ ...
    'Sample-size dependence of the designated primary ' ...
    'analysis only: LSD-PCB difference connectome, covariate-adjusted, ' ...
    'positive network, across all six experience dimensions'];
record.entry_script = entryScript;
record.entry_script_sha256 = file_sha256(entryScript);
record.plan_function = fullfile(fileparts(entryScript), ...
    'learning_curve_plan.m');
record.summary_builder = fullfile(fileparts(entryScript), ...
    'summarize_learning_curves.m');
record.established_engine = engineFile;
record.established_engine_sha256 = file_sha256(engineFile);
record.engine_provenance = ...
    'Publication copy of the validated optimized sample-size CPM engine.';
record.input_path = inputPath;
record.output_path = outputPath;
record.results_text = rawCsv;
record.results_summary = summaryCsv;
record.results_mat = checkpointMat;
record.run_log = runLog;
record.sample_plan = planPath;
record.behaviors = outcomes;
record.behavior_order = 'Table 1 order (GDE, MEQ30, VRS, OBN, BDE, AED)';
record.behavior_test_sets = 'within-LSD held-out folds only';
record.connectome_forms = {connectomeForm};
record.covariate_setting = covariateSetting;
record.analysis_grid = [ ...
    'six outcomes x LSD-PCB difference connectome x covariate-adjusted ' ...
    'only x outcome-appropriate sample sizes x 100 repetitions'];
record.standard_sample_sizes = standardSizes;
record.meq30_sample_sizes = meqSizes;
record.repetitions_per_sample_size = repetitions;
record.sampling_scheme = [ ...
    'nested study-stratified random subsampling'];
record.sample_rng = [ ...
    'standard rng(810000+repetition,twister); MEQ30 ' ...
    'rng(820000+repetition,twister)'];
record.cv_rng = sprintf('rng(cvSeed,%s), cvSeed=%d:%d', ...
    plan.cv_generator, plan.cv_seeds(1), plan.cv_seeds(end));
record.subject_subset = sprintf([ ...
    'GDE/VRS/OBN/BDE/AED use LSD rows 1:67; MEQ30 uses the established ' ...
    'subset rows %d:%d (maximum N=%d)'], meqSubset(1), meqSubset(end), ...
    numel(meqSubset));
record.difference_behavior_fallback = [ ...
    'use *_difference when present; otherwise use the standard behavior ' ...
    'vector, which is arithmetically identical because all placebo ' ...
    'ratings were zero (BDE and AED)'];
record.covariate_file = fullfile(inputPath, 'covars_indicator.mat');
record.covariate_columns = '[study_2, study_3, sex, age, mean_FD]';
record.covariate_coding_note = [ ...
    'Indicator-coded study with study 1 as reference; the raw numeric ' ...
    'study column in covars.mat is used only to build stratified samples'];
record.meq30_covariate_note = [ ...
    'The rows 20:67 subset contains only studies 2 and 3, so the study_2 ' ...
    'indicator is dropped for MEQ30 and study 2 becomes the reference ' ...
    'category. MEQ30 therefore uses 4 covariate columns ' ...
    '[study_3, sex, age, mean_FD] and a full-rank design. All other ' ...
    'outcomes use all 5 columns.'];
record.meq30_covariate_columns = '[study_3, sex, age, mean_FD]';
record.no_iterations = noIterations;
record.permutation_scheme = sprintf([ ...
    'observed plus shuffled behaviors; RandStream %s seed ' ...
    'iteration+%d'], rngSettings.internal_permutation_generator, ...
    rngSettings.permutation_seed_offset);
record.edge_threshold = edgeThreshold;
record.internal_folds = kFolds;
record.minimum_sample_size_note = [ ...
    'Configured sample sizes are recorded in standard_sample_sizes and meq30_sample_sizes'];
record.correlation_types = ...
    'GDE/BDE/AED Spearman; OBN/VRS/MEQ30 Pearson';
record.performance_correlation = 'Pearson observed-versus-predicted r';
record.summary_scope = [ ...
    'positive CPM network is the reported model: median r, 2.5th/97.5th ' ...
    'repetition percentiles, median raw p, and raw-p<0.05 frequency; the ' ...
    'raw file and summary retain combined and negative networks too'];
record.multiple_comparisons_correction = ...
    'none; permutation p-values are raw and detection frequency is descriptive';
record.parallel_workers_requested = workers;
record.expected_result_rows = expectedRows;
record.completed_result_rows = 0;
record.matlab_release = version('-release');
record.matlab_version = version;
end


function value = timestamp_now
value = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS'));
end
