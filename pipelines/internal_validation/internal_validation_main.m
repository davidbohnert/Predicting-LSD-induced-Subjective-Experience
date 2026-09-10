function results = internal_validation_main(cfg, presetName, options)
%INTERNAL_VALIDATION_MAIN Publication-facing internal CPM validation.
%
% RESULTS = INTERNAL_VALIDATION_MAIN(CFG, PRESETNAME, OPTIONS) runs a named
% preset from paper_analysis_config. This consolidates the primary and
% supplementary internal-validation specifications without outcome-specific
% launch scripts. OPTIONS may override no_iterations, workers, output_path,
% and use_parallel.

arguments
    cfg (1,1) struct
    presetName (1,1) string = "primary"
    options.no_iterations (1,1) double = NaN
    options.workers (1,1) double = NaN
    options.output_path (1,1) string = ""
    options.use_parallel (1,1) logical = true
    options.bootstrap_iterations (1,1) double = NaN
    options.connectome_directory (1,1) string = ""
end

pipelineDir = fileparts(mfilename('fullpath'));
repositoryRoot = fileparts(fileparts(pipelineDir));
addpath(pipelineDir, fullfile(repositoryRoot, 'scripts', 'common'), ...
    fullfile(repositoryRoot, 'scripts', 'utilities'));

assert(isfield(cfg, 'presets') && isfield(cfg.presets, presetName), ...
    'CPM:Publication:UnknownPreset', 'Unknown internal preset: %s', presetName);
preset = cfg.presets.(presetName);
if isnan(options.no_iterations), noIterations = cfg.runtime.main_iterations;
else, noIterations = options.no_iterations; end
if isnan(options.workers), workers = cfg.runtime.workers;
else, workers = options.workers; end
validateattributes(noIterations, {'numeric'}, {'scalar','integer','>=',2});
validateattributes(workers, {'numeric'}, {'scalar','integer','>=',0});
if options.output_path == ""
    outputPath = fullfile(cfg.output_root, 'internal_validation', presetName);
else
    outputPath = char(options.output_path);
end
if ~isfolder(outputPath), mkdir(outputPath); end
outputText = fullfile(outputPath, 'internal_validation_results.txt');
if options.connectome_directory == ""
    connectomeDirectory = fullfile(cfg.input_root, 'connectomes');
else
    connectomeDirectory = char(options.connectome_directory);
end

inputFiles = cellfun(@(stem) fullfile(connectomeDirectory, [stem '.mat']), ...
    preset.connectomes, 'UniformOutput',false);
[inputFiles, codeFiles] = cpm_analysis_dependencies(cfg, {pipelineDir}, inputFiles);
settings = struct('preset',preset, 'outcomes',cfg.outcomes, ...
    'model',rmfield(cfg.model, {'covariate_file','grouping_file'}), ...
    'rng',cfg.rng, 'iterations',noIterations);
cpm_checkpoint_guard(outputPath, settings, inputFiles, codeFiles);
checkpointPath = fullfile(outputPath, 'internal_validation_checkpoint.mat');
started = tic;
record = struct('status','RUNNING', 'analysis','internal_validation', ...
    'preset',char(presetName), 'entry_script',[mfilename('fullpath') '.m'], ...
    'input_path',cfg.input_root, 'output_path',outputPath, ...
    'no_iterations',noIterations, 'parallel_workers_requested',workers, ...
    'covariates',['full-rank indicator coding selected for represented ' ...
        'studies; see MAT details for each analysis'], ...
    'randomization_scheme', sprintf([ ...
        'folds rng(%d,%s); permutations iteration+%d,%s'], ...
        cfg.rng.cv_seed, cfg.rng.cv_generator, ...
        cfg.rng.permutation_seed_offset, ...
        cfg.rng.internal_permutation_generator), ...
    'completed_result_rows',0);
record.invocation_started = char(datetime('now'));
record.covariates_input_file = cfg.model.covariate_file;
record.grouping_file = cfg.model.grouping_file;
record.behaviors = preset.outcomes;
record.connectomes = preset.connectomes;
record.thresholds = preset.thresholds;
record.folds = preset.folds;
record.subject_subsets = 'configured outcome rows; MEQ30 20:67';
record.difference_behavior_fallback = 'prefer difference vector; standard equals difference when placebo is zero';
record.canonical_result_files = {checkpointPath, outputText, ...
    fullfile(outputPath,'internal_validation_results.csv'), ...
    fullfile(outputPath,'internal_validation_results.mat')};
record.parallel_workers_active = workers * double(options.use_parallel);
runRecordPath = save_analysis_run_record(outputPath, '', record);

if options.use_parallel && workers > 0
    pool = gcp('nocreate');
    if ~isempty(pool) && pool.NumWorkers ~= workers, delete(pool); pool = []; end
    if isempty(pool), parpool('Processes', workers); end
end

adjustments = logical(preset.adjustment);
thresholds = preset.thresholds;
foldValues = preset.folds;
outcomes = select_outcomes(cfg.outcomes, preset.outcomes);
if isfield(preset, 'subject_rows'), forcedRows = preset.subject_rows;
else, forcedRows = []; end
invertCorrelation = isfield(preset, 'invert_correlation_type') && ...
    logical(preset.invert_correlation_type);
useStudyLoso = isfield(preset, 'fold_source') && ...
    strcmp(preset.fold_source, 'study_labels');

covarsFull = load_single_numeric(cfg.model.covariate_file);
assert(isequal(size(covarsFull), [67 5]), ...
    'CPM:Publication:Covariates', ...
    'covars_indicator.mat must be 67-by-5.');
grouping = load_single_numeric(cfg.model.grouping_file);
assert(isequal(size(grouping), [67 4]), ...
    'CPM:Publication:Grouping', 'covars.mat must be 67-by-4.');
studyLabelsFull = grouping(:, 1);

resultRows = cell(0, 19);
details = struct([]);
completedUnits = {};
if isfile(checkpointPath)
    saved = load(checkpointPath);
    resultRows = saved.resultRows;
    details = saved.details;
    completedUnits = saved.completedUnits;
end
detailIndex = numel(details);
for connectomeIndex = 1:numel(preset.connectomes)
    connectomeStem = preset.connectomes{connectomeIndex};
    matsFull = load_single_numeric(fullfile(connectomeDirectory, ...
        [connectomeStem '.mat']));
    validate_connectome(matsFull, connectomeStem);

    for outcomeIndex = 1:numel(outcomes)
        outcome = outcomes(outcomeIndex);
        if isempty(forcedRows), subjectRows = outcome.subject_rows;
        else, subjectRows = forcedRows; end
        mats = matsFull(:, :, subjectRows);
        behaviorStem = behavior_for_connectome(outcome, connectomeStem, ...
            cfg.input_root);
        behavior = load_single_numeric(fullfile(cfg.input_root, 'behav', ...
            [behaviorStem '.mat']));
        behavior = behavior(:);
        if numel(behavior) == 67 && numel(subjectRows) < 67
            behavior = behavior(subjectRows);
        end
        assert(numel(behavior) == numel(subjectRows), ...
            'CPM:Publication:BehaviorRows', ...
            '%s does not match its configured subject subset.', behaviorStem);

        corrType = outcome.correlation_type;
        if invertCorrelation
            if strcmp(corrType, 'Pearson'), corrType = 'Spearman';
            else, corrType = 'Pearson'; end
        end

        for adjustment = adjustments
            if adjustment
                if useStudyLoso
                    % Study-LOSO chooses the reference category separately
                    % inside each training fold.
                    covars = covarsFull(subjectRows, :);
                    covariateNames = {'fold-specific full-rank coding'};
                    referenceStudy = NaN;
                else
                    [covars, covariateNames, referenceStudy] = ...
                        cpm_select_covariates(covarsFull(subjectRows, :), ...
                            studyLabelsFull(subjectRows));
                end
            else
                covars = [];
                covariateNames = {};
                referenceStudy = NaN;
            end
            for foldValue = foldValues
                unit = sprintf('%s/%s/%d/%g', connectomeStem, outcome.label, adjustment, foldValue);
                if ismember(unit, completedUnits), continue; end
                if isinf(foldValue), k = numel(subjectRows); else, k = foldValue; end
                if useStudyLoso
                    assert(numel(thresholds) == 1 && adjustment, ...
                        'CPM:Publication:StudyLOSOConfiguration', ...
                        'Study LOSO requires one threshold and adjustment.');
                    [studyResult, runDetails] = ...
                        permutation_test_cv_study_loso_positive( ...
                            connectomeStem, behaviorStem, mats, behavior, ...
                            studyLabelsFull(subjectRows), covars, thresholds, ...
                            noIterations, corrType, ...
                            workers * double(options.use_parallel), cfg.rng);
                    analysisLabel = sprintf('%s vs %s - with covariates', ...
                        connectomeStem, behaviorStem);
                    rows = {analysisLabel, numel(unique( ...
                        studyLabelsFull(subjectRows))), thresholds, ...
                        NaN, NaN, studyResult.r, studyResult.p, NaN, NaN, ...
                        NaN, NaN, studyResult.mse, NaN, NaN, NaN};
                    consensusResults = {};
                else
                    if options.use_parallel, engineWorkers = workers;
                    else, engineWorkers = 0; end
                    [rows, runDetails] = permutation_test_cv_optimized( ...
                        connectomeStem, behaviorStem, mats, behavior, k, ...
                        thresholds, covars, outputText, {}, noIterations, ...
                        corrType, engineWorkers, cfg.rng);
                    rng(cfg.rng.cv_seed, cfg.rng.cv_generator);
                    fixedPartition = cvpartition(numel(subjectRows), ...
                        'KFold', k);
                    consensusResults = cell(numel(thresholds), 1);
                    for thresholdIndex = 1:numel(thresholds)
                        consensusResults{thresholdIndex} = ...
                            cpm_build_consensus_masks(mats, behavior, covars, ...
                                corrType, fixedPartition, ...
                                thresholds(thresholdIndex), ...
                                cfg.model.consensus_fraction, engineWorkers);
                    end
                end
                for rowIndex = 1:size(rows, 1)
                    rows{rowIndex, 16} = char(presetName);
                    rows{rowIndex, 17} = noIterations;
                    if isempty(consensusResults)
                        rows{rowIndex, 18} = NaN;
                        rows{rowIndex, 19} = NaN;
                    else
                        rows{rowIndex, 18} = nnz(triu( ...
                            consensusResults{rowIndex}.observed_pos_mask, 1));
                        rows{rowIndex, 19} = nnz(triu( ...
                            consensusResults{rowIndex}.observed_neg_mask, 1));
                    end
                end
                resultRows = [resultRows; rows]; %#ok<AGROW>
                detailIndex = detailIndex + 1;
                details(detailIndex).outcome = outcome.label; %#ok<AGROW>
                details(detailIndex).connectome = connectomeStem;
                details(detailIndex).adjusted = adjustment;
                details(detailIndex).subject_rows = subjectRows;
                details(detailIndex).correlation_type = corrType;
                details(detailIndex).covariate_columns = covariateNames;
                details(detailIndex).study_reference = referenceStudy;
                if useStudyLoso
                    details(detailIndex).covariate_design_rank = NaN;
                else
                    details(detailIndex).covariate_design_rank = ...
                        rank([ones(numel(subjectRows), 1), covars]);
                end
                details(detailIndex).engine_details = runDetails;
                details(detailIndex).consensus = consensusResults;
                if strcmp(presetName, "primary") && ~isempty(consensusResults)
                    maskDirectory = fullfile(outputPath, 'consensus_masks');
                    if ~isfolder(maskDirectory), mkdir(maskDirectory); end
                    positiveMask = consensusResults{1}.observed_pos_mask;
                    negativeMask = consensusResults{1}.observed_neg_mask;
                    maskMetadata = struct('outcome',outcome.label, ...
                        'connectome',connectomeStem, 'threshold',thresholds(1), ...
                        'folds',k, 'consensus_fraction', ...
                        cfg.model.consensus_fraction, ...
                        'correlation_type',corrType, ...
                        'subject_rows',subjectRows, ...
                        'covariate_columns',{covariateNames}, ...
                        'study_reference',referenceStudy);
                    save(fullfile(maskDirectory, ...
                        ['positive_consensus_' outcome.label '.mat']), ...
                        'positiveMask','negativeMask','maskMetadata','-v7');
                end
                completedUnits{end+1} = unit; %#ok<AGROW>
                payload = struct('resultRows',{resultRows}, 'details',details, ...
                    'completedUnits',{completedUnits});
                cpm_atomic_save(checkpointPath, payload);
                record.completed_result_rows = size(resultRows,1);
                record.covariate_columns = covariateNames;
                record.study_reference = referenceStudy;
                save_analysis_run_record(outputPath, runRecordPath, record);
            end
        end
    end
end

header = {'Analysis','k','threshold','r_combined','p_combined', ...
    'r_positive','p_positive','r_negative','p_negative','mse_combined', ...
    'mse_p_combined','mse_positive','mse_p_positive','mse_negative', ...
    'mse_p_negative','preset','iterations','positive_consensus_edges', ...
    'negative_consensus_edges'};
results = cell2table(resultRows, 'VariableNames', header);
writetable(results, fullfile(outputPath, 'internal_validation_results.csv'));
if strcmp(presetName, "primary")
    if isnan(options.bootstrap_iterations)
        bootstrapIterations = noIterations;
    else
        bootstrapIterations = options.bootstrap_iterations;
    end
    validateattributes(bootstrapIterations, {'numeric'}, ...
        {'scalar','integer','>=',2});
    primaryTable = build_primary_internal_table( ...
        results, details, bootstrapIterations);
    writetable(primaryTable, fullfile(outputPath, ...
        'table_1_internal_validation.csv'));
else
    primaryTable = table;
end
save(fullfile(outputPath, 'internal_validation_results.mat'), ...
    'results', 'details', 'primaryTable', 'cfg', 'presetName', '-v7.3');

record.status = 'COMPLETE';
record.runtime_seconds = toc(started);
record.invocation_finished = char(datetime('now'));
record.covariate_columns = arrayfun(@(d) ...
    [d.outcome ': ' strjoin(d.covariate_columns, ', ')], details, ...
    'UniformOutput',false);
record.completed_result_rows = height(results);
if strcmp(presetName, "primary")
    record.bootstrap_iterations = bootstrapIterations;
    record.bootstrap_rng = ...
        'rng(1300000 + sample size,twister); fixed OOF predictions';
end
save_analysis_run_record(outputPath, runRecordPath, record);
end

function selected = select_outcomes(allOutcomes, labels)
selected = allOutcomes(ismember({allOutcomes.label}, labels));
assert(numel(selected) == numel(labels), ...
    'CPM:Publication:OutcomeSelection', 'A configured outcome is missing.');
end

function stem = behavior_for_connectome(outcome, connectomeStem, inputRoot)
if contains(lower(connectomeStem), 'difference')
    preferred = outcome.difference_behavior;
    if isfile(fullfile(inputRoot, 'behav', [preferred '.mat']))
        stem = preferred;
    else
        stem = outcome.standard_behavior;
    end
else
    stem = outcome.standard_behavior;
end
end

function value = load_single_numeric(pathName)
assert(isfile(pathName), 'CPM:Publication:MissingInput', ...
    'Missing input: %s', pathName);
loaded = load(pathName);
names = fieldnames(loaded);
assert(isscalar(names) && isnumeric(loaded.(names{1})), ...
    'CPM:Publication:InputFormat', ...
    '%s must contain one numeric variable.', pathName);
value = loaded.(names{1});
assert(all(isfinite(value), 'all'), 'CPM:Publication:NonfiniteInput', ...
    '%s contains NaN or Inf.', pathName);
end

function validate_connectome(mats, name)
assert(ndims(mats) == 3 && size(mats,1) == size(mats,2) && ...
    size(mats,3) == 67, 'CPM:Publication:ConnectomeShape', ...
    '%s must be node-by-node-by-67.', name);
end
