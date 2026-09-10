function results = cross_drug_loso_main(cfg, options)
%CROSS_DRUG_LOSO_MAIN Subject-wise cross-drug CPM with threshold sensitivity.
% Train on LSD difference connectomes after excluding the matching participant.
% Pool the two other-drug sessions within each pharmacological class.
% Completed outcome/class cells are checkpointed for compatible resumes.

arguments
    cfg (1,1) struct
    options.no_iterations (1,1) double = NaN
    options.workers (1,1) double = NaN
    options.output_path (1,1) string = ""
    options.scales cell = {'GDE','MEQ30','VRS'}
    options.thresholds (1,:) double = [0.01 0.05 0.005]
    options.connectome_form (1,1) string = "difference"
    options.connectome_directory (1,1) string = ""
    options.connectome_suffix (1,1) string = ""
end


scriptRoot = fileparts(mfilename('fullpath'));
repositoryRoot = fileparts(fileparts(scriptRoot));
inputPath = cfg.input_root;
if options.connectome_directory == ""
    connectomeDirectory = fullfile(inputPath, 'connectomes');
else
    connectomeDirectory = char(options.connectome_directory);
end
crosswalkFile = cfg.crosswalk_file;
utilitiesPath = fullfile(repositoryRoot, 'scripts', 'utilities');
commonPath = fullfile(repositoryRoot, 'scripts', 'common');

if options.output_path == ""
    outputPath = fullfile(cfg.output_root, 'cross_drug_loso');
else
    outputPath = char(options.output_path);
end
auditPath = fullfile(outputPath, 'audit');

addpath(scriptRoot, utilitiesPath, commonPath);
if ~isfolder(outputPath), mkdir(outputPath); end


runTag = 'cross_drug_loso';
outputMat = fullfile(outputPath, [runTag '_checkpoint.mat']);
outputCsv = fullfile(outputPath, [runTag '_results.csv']);
outputXlsx = fullfile(outputPath, [runTag '_results.xlsx']);
outputTxt = fullfile(outputPath, [runTag '_results.txt']);
runLog = fullfile(outputPath, [runTag '_run_log.txt']);

SCALES = options.scales;
CONNECTOME_FORM = char(validatestring(options.connectome_form, ...
    ["standard","difference"]));
COVARIATE_SETTING = 'with_covariates';
MAIN_THRESHOLD = 0.01;
P_THRESHOLDS = options.thresholds;
if isnan(options.no_iterations), NO_ITERATIONS = cfg.runtime.main_iterations;
else, NO_ITERATIONS = options.no_iterations; end
if isnan(options.workers), POOL_WORKERS = cfg.runtime.workers;
else, POOL_WORKERS = options.workers; end
MEQ30_SUBSET = (20:67)';
SUBJECT_COUNT = 67;

cohortConfigs = [ ...
    struct('cohort', "LAM", 'combinedStem', "all_amphs", ...
           'drugAStem', "amphetamine", 'drugACrosswalk', "amphetamine", ...
           'drugBStem', "MDMA", 'drugBCrosswalk', "mdma", ...
           'expectedSubjects', 23), ...
    struct('cohort', "LPM", 'combinedStem', "other_psych", ...
           'drugAStem', "psilocybin", 'drugACrosswalk', "psilocybin", ...
           'drugBStem', "mescaline", 'drugBCrosswalk', "mescaline", ...
           'expectedSubjects', 25) ...
];

header = result_header();

connectomeFiles = dir(fullfile(connectomeDirectory, ['*' char(options.connectome_suffix) '.mat']));
extraInputs = arrayfun(@(f) fullfile(f.folder,f.name),connectomeFiles,'UniformOutput',false);
[inputs, code] = cpm_analysis_dependencies(cfg, {scriptRoot}, [extraInputs' {crosswalkFile}]);
settings = struct('scales',{SCALES}, 'thresholds',P_THRESHOLDS, ...
    'iterations',NO_ITERATIONS, 'form',CONNECTOME_FORM, ...
    'suffix',char(options.connectome_suffix), 'rng',cfg.rng, 'outcomes',cfg.outcomes);
cpm_checkpoint_guard(outputPath, settings, inputs, code);
if ~isfolder(auditPath), mkdir(auditPath); end
diary(runLog);
diary on;
diaryCleanup = onCleanup(@() diary('off'));

fprintf('Publication cross-drug LOSO analysis\n');
fprintf('Invocation started: %s\n', timestamp_now());
fprintf('Scales: %s\n', strjoin(SCALES, ', '));
fprintf('Connectome form: %s; covariates: %s\n', ...
    CONNECTOME_FORM, COVARIATE_SETTING);
fprintf('Thresholds: %s (main %g)\n', mat2str(P_THRESHOLDS), MAIN_THRESHOLD);
fprintf('Iterations: %d; requested workers: %d\n', NO_ITERATIONS, POOL_WORKERS);
fprintf('Output: %s\n', outputPath);

entryScript = [mfilename('fullpath') '.m'];
runRecord = struct;
runRecord.status = 'RUNNING';
runRecord.invocation_started = timestamp_now();
runRecord.analysis = 'cross_drug_loso';
runRecord.purpose = ['Publication cross-drug edge-retention threshold ' ...
    'analysis for the primary covariate-adjusted LSD-PCB ' ...
    'difference model, subject-wise LOSO'];
runRecord.entry_script = entryScript;
runRecord.permutation_engine = fullfile(scriptRoot, ...
    'permutation_test_cross_drug_loso.m');
runRecord.input_path = inputPath;
runRecord.crosswalk = crosswalkFile;
runRecord.output_path = outputPath;
runRecord.results_mat = outputMat;
runRecord.results_csv = outputCsv;
runRecord.results_workbook = outputXlsx;
runRecord.results_text = outputTxt;
runRecord.run_log = runLog;
runRecord.audit_path = auditPath;
runRecord.behaviors = SCALES;
runRecord.test_sets = ['psilocybin,mescaline,other_psych (LPM); ' ...
    'amphetamine,MDMA,all_amphs (LAM)'];
runRecord.connectome_forms = CONNECTOME_FORM;
runRecord.covariate_settings = ...
    'covars_indicator.mat [study_2, study_3, sex, age, mean_FD]';
runRecord.main_edge_threshold = MAIN_THRESHOLD;
runRecord.control_edge_thresholds = mat2str(setdiff(P_THRESHOLDS, MAIN_THRESHOLD));
runRecord.consensus_thresholds = 'not applicable to direct LOSO';
runRecord.no_iterations = NO_ITERATIONS;
runRecord.feature_selection_correlations = ...
    'Pearson for OBN/VRS/MEQ30, Spearman for BDE/GDE/AED';
runRecord.performance_correlation = 'Pearson';
runRecord.feature_selection = ...
    'direct on all n-1 eligible LSD training subjects within each fold';
runRecord.subject_subset = ...
    'LSD rows 1:67; MEQ30 uses the established 20:67 subset';
runRecord.holdout_rule = ['remove the matching LSD subject before every ' ...
    'supervised operation; predict both other-drug sessions'];
runRecord.permutation_scheme = ['within-cohort whole-subject behavior ' ...
    'profiles: matching LSD and both test-drug outcomes move together'];
runRecord.random_seeds = sprintf( ...
    'LOSO subject-profile permutations rng(iteration+%d,%s)', ...
    cfg.rng.permutation_seed_offset, cfg.rng.loso_permutation_generator);
runRecord.multiple_comparisons_correction = [ ...
    'raw p-values retained in the complete output; Table 2 applies BH ' ...
    'across GDE, MEQ30, and VRS separately within each drug class'];
runRecord.difference_behavior_fallback = 'true';
runRecord.parallel_workers_requested = POOL_WORKERS;
runRecordPath = save_analysis_run_record(outputPath, '', runRecord);

[covarsAll, covarVariable] = load_first_var_with_name( ...
    cfg.model.covariate_file);
assert(isnumeric(covarsAll) && isequal(size(covarsAll), [SUBJECT_COUNT, 5]) ...
    && all(isfinite(covarsAll(:))), ...
    'CPM:SuppLOSO:InvalidCovariates', ...
    'covars_indicator.mat must contain a finite 67-by-5 matrix.');
assert(all(ismember(covarsAll(:, 1:3), [0 1]), 'all'), ...
    'CPM:SuppLOSO:CovariateCodingMismatch', ...
    ['Columns 1-3 must be 0/1 indicators (study_2, study_3, sex). ' ...
     'Numeric study labels mean the legacy covars.mat was loaded.']);
fprintf('Covariate variable: %s (%d-by-%d)\n', ...
    covarVariable, size(covarsAll, 1), size(covarsAll, 2));
[groupingAll, ~] = load_first_var_with_name(cfg.model.grouping_file);
assert(isequal(size(groupingAll), [SUBJECT_COUNT, 4]), ...
    'CPM:SuppLOSO:InvalidGrouping', ...
    'covars.mat must contain a 67-by-4 grouping matrix.');
studyLabelsAll = groupingAll(:, 1);

crosswalk = readtable(crosswalkFile, 'TextType', 'string');
validate_crosswalk_schema(crosswalk);

if strcmp(CONNECTOME_FORM, 'difference')
    trainMatrixStem = 'LSD_difference';
    matrixSuffix = '_difference';
else
    trainMatrixStem = 'LSD';
    matrixSuffix = '';
end
trainMatrixFile = [trainMatrixStem char(options.connectome_suffix) '.mat'];
matsTrainFull = load_first_var(fullfile(connectomeDirectory, trainMatrixFile));
validate_connectome(matsTrainFull, trainMatrixStem, SUBJECT_COUNT);

cells = build_cell_plan(SCALES, cohortConfigs, inputPath, CONNECTOME_FORM);
expectedRows = 3 * numel(P_THRESHOLDS) * numel(cells);
runRecord.expected_result_rows = expectedRows;
fprintf('Planned cells: %d; expected result rows: %d\n', ...
    numel(cells), expectedRows);
report_skipped_cells(SCALES, cohortConfigs, inputPath);

parameters = parameter_sheet(NO_ITERATIONS, POOL_WORKERS, P_THRESHOLDS, ...
    MAIN_THRESHOLD, expectedRows, SCALES, cfg.rng);

results = load_results_checkpoint(outputMat, header);
runRecord.completed_result_rows = max(size(results, 1) - 1, 0);
runRecordPath = save_analysis_run_record( ...
    outputPath, runRecordPath, runRecord);

pool = gcp('nocreate');
if POOL_WORKERS > 0
    if ~isempty(pool) && pool.NumWorkers ~= POOL_WORKERS
        delete(pool);
        pool = [];
    end
    if isempty(pool), pool = parpool('Processes', POOL_WORKERS); end
    assert(pool.NumWorkers == POOL_WORKERS, ...
        'CPM:SuppLOSO:WorkerCount', ...
        'Expected %d workers; started %d.', POOL_WORKERS, pool.NumWorkers);
    activeWorkers = pool.NumWorkers;
else
    activeWorkers = 0;
end
runRecord.parallel_workers_active = activeWorkers;
runRecordPath = save_analysis_run_record( ...
    outputPath, runRecordPath, runRecord);
fprintf('Parallel workers active: %d\n', activeWorkers);

analysisStart = tic;
for cellIdx = 1:numel(cells)
    spec = cells(cellIdx);
    config = spec.config;
    analysisId = spec.analysisId;
    auditFile = fullfile(auditPath, [analysisId '.mat']);

    if analysis_is_complete(results, analysisId, P_THRESHOLDS, config) ...
            && isfile(auditFile)
        fprintf('SKIP completed: %s\n', analysisId);
        continue
    end

    outcomeSpec = cfg.outcomes(strcmp({cfg.outcomes.label},spec.scale));
    assert(isscalar(outcomeSpec), 'CPM:Config:Outcome', 'Unknown outcome.');
    corrType = outcomeSpec.correlation_type;

    if strcmp(spec.scale, 'MEQ30')
        trainGlobalRows = MEQ30_SUBSET;
    else
        trainGlobalRows = (1:SUBJECT_COUNT)';
    end
    matsTrain = matsTrainFull(:, :, trainGlobalRows);
    [covarsTrain, covariateNames, referenceStudy] = ...
        cpm_select_covariates(covarsAll(trainGlobalRows, :), ...
            studyLabelsAll(trainGlobalRows));
    nTrainPerFold = size(matsTrain, 3) - 1;

    standardTrainBehavior = ['LSD_' spec.scale];
    if strcmp(CONNECTOME_FORM, 'difference')
        [behavTrain, trainBehaviorSource] = load_behavior_with_fallback( ...
            inputPath, [standardTrainBehavior '_difference'], ...
            standardTrainBehavior);
    else
        behavTrain = load_behavior(fullfile(inputPath,'behav', ...
            [standardTrainBehavior '.mat']), standardTrainBehavior);
        trainBehaviorSource = standardTrainBehavior;
    end
    assert(numel(behavTrain) == numel(trainGlobalRows), ...
        'CPM:SuppLOSO:TrainBehaviorCount', ...
        '%s must contain %d values.', ...
        trainBehaviorSource, numel(trainGlobalRows));
    trainCohorts = cohort_labels_from_lsd_rows(trainGlobalRows);

    [heldoutGlobalRows, anonymousSubjectIndex] = ...
        validate_cohort_crosswalk(crosswalk, config);
    [foundHeldout, heldoutTrainIdx] = ismember( ...
        heldoutGlobalRows, trainGlobalRows);
    assert(all(foundHeldout), ...
        'CPM:SuppLOSO:HeldoutMissingFromTraining', ...
        'A target-cohort subject is absent from LSD training.');
    heldoutTrainIdx = heldoutTrainIdx(:);

    drugAMatrixStem = char(config.drugAStem + matrixSuffix);
    drugBMatrixStem = char(config.drugBStem + matrixSuffix);
    matsTestA = load_first_var(fullfile(connectomeDirectory, ...
        [drugAMatrixStem char(options.connectome_suffix) '.mat']));
    matsTestB = load_first_var(fullfile(connectomeDirectory, ...
        [drugBMatrixStem char(options.connectome_suffix) '.mat']));
    validate_connectome(matsTestA, drugAMatrixStem, config.expectedSubjects);
    validate_connectome(matsTestB, drugBMatrixStem, config.expectedSubjects);
    validate_single_drug_connectomes(connectomeDirectory, config, ...
        [matrixSuffix char(options.connectome_suffix)], matsTestA, matsTestB);

    standardTestA = char(config.drugAStem + "_" + spec.scale);
    standardTestB = char(config.drugBStem + "_" + spec.scale);
    if strcmp(CONNECTOME_FORM, 'difference')
        [behavTestA, testBehaviorSourceA] = load_behavior_with_fallback( ...
            inputPath, [standardTestA '_difference'], standardTestA);
        [behavTestB, testBehaviorSourceB] = load_behavior_with_fallback( ...
            inputPath, [standardTestB '_difference'], standardTestB);
    else
        behavTestA = load_behavior(fullfile(inputPath,'behav', ...
            [standardTestA '.mat']), standardTestA);
        behavTestB = load_behavior(fullfile(inputPath,'behav', ...
            [standardTestB '.mat']), standardTestB);
        testBehaviorSourceA = standardTestA;
        testBehaviorSourceB = standardTestB;
    end
    [testBehaviorSourceCombined, combinedChecked] = ...
        validate_single_drug_behaviors(inputPath, config, spec.scale, ...
            behavTestA, behavTestB, CONNECTOME_FORM);

    fprintf('\n========================================\n');
    fprintf('[%d/%d] START %s at %s\n', ...
        cellIdx, numel(cells), analysisId, timestamp_now());
    fprintf('Training: %s rows %d:%d (%s), n/fold=%d\n', ...
        trainMatrixStem, trainGlobalRows(1), trainGlobalRows(end), ...
        trainBehaviorSource, nTrainPerFold);
    fprintf('Test behavior: %s / %s / %s\n', ...
        testBehaviorSourceA, testBehaviorSourceB, testBehaviorSourceCombined);
    fprintf('Correlation: %s; thresholds: %s\n', ...
        corrType, mat2str(P_THRESHOLDS));
    fprintf('========================================\n');

    started = tic;
    output = permutation_test_cross_drug_loso( ...
        matsTrain, behavTrain, covarsTrain, trainGlobalRows, trainCohorts, ...
        heldoutTrainIdx, config.cohort, matsTestA, behavTestA, ...
        matsTestB, behavTestB, P_THRESHOLDS, corrType, ...
        NO_ITERATIONS, POOL_WORKERS > 0, cfg.rng);
    runtimeSeconds = toc(started);

    assert(output.audit.n_train_per_fold == nTrainPerFold, ...
        'CPM:SuppLOSO:TrainPerFold', ...
        'Expected %d training subjects per fold; engine used %d.', ...
        nTrainPerFold, output.audit.n_train_per_fold);

    audit = output.audit;
    audit.analysis_id = analysisId;
    audit.connectome_form = CONNECTOME_FORM;
    audit.behavior_scale = spec.scale;
    audit.training_matrix = trainMatrixStem;
    audit.training_behavior_source = trainBehaviorSource;
    audit.subject_subset = trainGlobalRows;
    audit.covariate_setting = COVARIATE_SETTING;
    audit.covariate_variable = covarVariable;
    audit.covariate_column_order = strjoin(covariateNames, ', ');
    audit.study_reference = referenceStudy;
    audit.covariate_design_rank = ...
        rank([ones(numel(trainGlobalRows), 1), covarsTrain]);
    audit.feature_selection_correlation = corrType;
    audit.drug_a = config.drugAStem;
    audit.drug_b = config.drugBStem;
    audit.combined_test_set = config.combinedStem;
    audit.drug_a_behavior_source = testBehaviorSourceA;
    audit.drug_b_behavior_source = testBehaviorSourceB;
    audit.combined_behavior_source = testBehaviorSourceCombined;
    audit.combined_behavior_cross_checked = combinedChecked;
    audit.anonymous_subject_index = anonymousSubjectIndex;
    audit.runtime_seconds = runtimeSeconds;
    save_audit_checkpoint(auditFile, audit);

    results = remove_analysis_rows(results, analysisId);
    for thresholdIdx = 1:numel(output.thresholds)
        threshold = output.thresholds(thresholdIdx);
        if threshold == MAIN_THRESHOLD
            specification = 'main';
        else
            specification = 'edge_threshold_control';
        end
        results = [results; ...
            build_result_row(header, analysisId, specification, ...
                trainMatrixStem, spec.scale, trainBehaviorSource, ...
                config.drugAStem, testBehaviorSourceA, config, ...
                CONNECTOME_FORM, COVARIATE_SETTING, corrType, threshold, ...
                NO_ITERATIONS, nTrainPerFold, config.expectedSubjects, ...
                true, output.drug_a(thresholdIdx)); ...
            build_result_row(header, analysisId, specification, ...
                trainMatrixStem, spec.scale, trainBehaviorSource, ...
                config.drugBStem, testBehaviorSourceB, config, ...
                CONNECTOME_FORM, COVARIATE_SETTING, corrType, threshold, ...
                NO_ITERATIONS, nTrainPerFold, config.expectedSubjects, ...
                true, output.drug_b(thresholdIdx)); ...
            build_result_row(header, analysisId, specification, ...
                trainMatrixStem, spec.scale, trainBehaviorSource, ...
                config.combinedStem, testBehaviorSourceCombined, config, ...
                CONNECTOME_FORM, COVARIATE_SETTING, corrType, threshold, ...
                NO_ITERATIONS, nTrainPerFold, ...
                2 * config.expectedSubjects, combinedChecked, ...
                output.combined(thresholdIdx))]; %#ok<AGROW>
    end

    checkpoint_results(results, parameters, outputMat, outputCsv, ...
        outputXlsx, outputTxt);
    runRecord.completed_result_rows = size(results, 1) - 1;
    runRecordPath = save_analysis_run_record( ...
        outputPath, runRecordPath, runRecord);
    fprintf('COMPLETE %s in %.1f seconds (%.3f h)\n', ...
        analysisId, runtimeSeconds, runtimeSeconds / 3600);
end

assert(size(results, 1) - 1 == expectedRows, ...
    'CPM:SuppLOSO:IncompleteGrid', ...
    'Expected %d result rows; found %d.', expectedRows, size(results, 1) - 1);
table2Outcomes = {'GDE','MEQ30','VRS'};
if all(ismember(table2Outcomes, SCALES)) && ismember(MAIN_THRESHOLD, P_THRESHOLDS)
    table2 = build_table_2(results, MAIN_THRESHOLD);
    writetable(table2, fullfile(outputPath, 'table_2_cross_drug.csv'));
else
    fprintf(['Table 2 requires threshold .01 and ' ...
        'all of GDE, MEQ30, and VRS.\n']);
end

runRecord.status = 'COMPLETE';
runRecord.invocation_finished = timestamp_now();
runRecord.runtime_seconds = toc(analysisStart);
runRecord.completed_result_rows = size(results, 1) - 1;
save_analysis_run_record(outputPath, runRecordPath, runRecord);

fprintf('\nCross-drug LOSO analysis COMPLETE at %s\n', ...
    timestamp_now());
fprintf('Completed result rows: %d\n', runRecord.completed_result_rows);
fprintf('Runtime: %.1f seconds (%.3f h)\n', ...
    runRecord.runtime_seconds, runRecord.runtime_seconds / 3600);
fprintf('Results workbook: %s\n', outputXlsx);
fprintf('Run record: %s\n', runRecordPath);
clear diaryCleanup
diary('off');
end


% ---------------------------------------------------------------- planning

function cells = build_cell_plan(scales, cohortConfigs, inputPath, connectomeForm)
%BUILD_CELL_PLAN One struct per engine call, skipping unmeasured cells.
cells = struct('analysisId', {}, 'scale', {}, 'config', {});
for cohortIdx = 1:numel(cohortConfigs)
    config = cohortConfigs(cohortIdx);
    for scaleIdx = 1:numel(scales)
        scale = scales{scaleIdx};
        if ~cohort_measures_scale(inputPath, config, scale)
            continue
        end
        cells(end + 1) = struct( ... %#ok<AGROW>
            'analysisId', sprintf('%s__%s__%s__with_covariates', ...
                connectomeForm, scale, config.cohort), ...
            'scale', scale, 'config', config);
    end
end
end


function tf = cohort_measures_scale(inputPath, config, scale)
%COHORT_MEASURES_SCALE True when both test drugs have this scale on file.
stemA = char(config.drugAStem + "_" + scale);
stemB = char(config.drugBStem + "_" + scale);
tf = behavior_exists(inputPath, stemA) && behavior_exists(inputPath, stemB);
end


function tf = behavior_exists(inputPath, stem)
tf = isfile(fullfile(inputPath, 'behav', [stem '.mat'])) || ...
    isfile(fullfile(inputPath, 'behav', [stem '_difference.mat']));
end


function report_skipped_cells(scales, cohortConfigs, inputPath)
for cohortIdx = 1:numel(cohortConfigs)
    config = cohortConfigs(cohortIdx);
    for scaleIdx = 1:numel(scales)
        scale = scales{scaleIdx};
        if ~cohort_measures_scale(inputPath, config, scale)
            fprintf(['NOT ASSESSED: %s was not measured in the %s cohort ' ...
                '(%s / %s); reported as not assessed, not as a null ' ...
                'result.\n'], scale, config.cohort, ...
                config.drugAStem, config.drugBStem);
        end
    end
end
end


% ----------------------------------------------------------------- results

function header = result_header
header = { ...
    'Analysis ID', 'Specification', 'Training Data', 'Behavior Scale', ...
    'Training Behavior Source', 'Testing Set', 'Testing Behavior Source', ...
    'Cohort', 'Connectome Form', 'Covariate Setting', ...
    'Feature Selection Correlation', 'Performance Correlation', ...
    'P Threshold', 'Iterations', 'N Train Per Fold', 'N Test', ...
    'Pooled Vector Cross Checked', ...
    'r_comb', 'p_comb', 'r_pos', 'p_pos', 'r_neg', 'p_neg', ...
    'mse_comb', 'mse_p_comb', 'mse_pos', 'mse_p_pos', ...
    'mse_neg', 'mse_p_neg', 'mse_baseline', 'P Value Note'};
end


function paperTable = build_table_2(results, mainThreshold)
outcomes = ["GDE"; "MEQ30"; "VRS"];
testSets = ["other_psych", "all_amphs"];
r = zeros(3,2); p = zeros(3,2);
body = results(2:end, :);
for outcomeIndex = 1:3
    for setIndex = 1:2
        matches = strcmp(body(:,4), outcomes(outcomeIndex)) & ...
            strcmp(body(:,6), testSets(setIndex)) & ...
            cellfun(@(x) isequal(x, mainThreshold), body(:,13));
        assert(nnz(matches) == 1, 'CPM:Reporting:CrossDrugMatch', ...
            'Expected one main row for %s in %s.', ...
            outcomes(outcomeIndex), testSets(setIndex));
        row = body(matches, :);
        r(outcomeIndex,setIndex) = row{20};
        p(outcomeIndex,setIndex) = row{21};
    end
end
q = [benjamini_hochberg(p(:,1)), benjamini_hochberg(p(:,2))];
paperTable = table(outcomes, r(:,1), p(:,1), q(:,1), ...
    r(:,2), p(:,2), q(:,2), 'VariableNames', ...
    {'Outcome','OtherPsychedelicsR','OtherPsychedelicsP', ...
     'OtherPsychedelicsQ','StimulantsR','StimulantsP','StimulantsQ'});
end


function row = build_result_row(header, analysisId, specification, ...
        trainMatrix, behaviorScale, trainBehaviorSource, testSet, ...
        testBehaviorSource, config, connectomeForm, covariateTag, ...
        corrType, pThreshold, noIterations, nTrainPerFold, nTest, ...
        crossChecked, result)
row = cell(1, numel(header));
row(:) = {NaN};
row{1} = analysisId;
row{2} = specification;
row{3} = trainMatrix;
row{4} = behaviorScale;
row{5} = trainBehaviorSource;
row{6} = char(testSet);
row{7} = testBehaviorSource;
row{8} = char(config.cohort);
row{9} = connectomeForm;
row{10} = covariateTag;
row{11} = corrType;
row{12} = 'Pearson';
row{13} = pThreshold;
row{14} = noIterations;
row{15} = nTrainPerFold;
row{16} = nTest;
if crossChecked
    row{17} = 'yes';
else
    row{17} = 'no pooled vector on file';
end
row{18} = result.r(1);
row{19} = result.p_r(1);
row{20} = result.r(2);
row{21} = result.p_r(2);
row{22} = result.r(3);
row{23} = result.p_r(3);
row{24} = result.mse(1);
row{25} = result.p_mse(1);
row{26} = result.mse(2);
row{27} = result.p_mse(2);
row{28} = result.mse(3);
row{29} = result.p_mse(3);
row{30} = result.baseline_mse;
row{31} = 'Raw permutation p-value';
end


function parameters = parameter_sheet(noIterations, workers, pThresholds, ...
        mainThreshold, expectedRows, scales, rngSettings)
parameters = { ...
    'Parameter', 'Value'; ...
    'Analysis', 'Cross-drug LOSO edge-threshold analysis'; ...
    'Behaviors', strjoin(scales, ', '); ...
    'Connectome form', 'LSD-PCB difference only (primary specification)'; ...
    'Covariates', ...
        'covars_indicator.mat [study_2, study_3, sex, age, mean_FD]'; ...
    'LSD subject subset', 'rows 1:67; MEQ30 uses rows 20:67'; ...
    'Main edge threshold', mainThreshold; ...
    'Thresholds', mat2str(pThresholds); ...
    'Consensus thresholds', ...
        'Not applicable under direct LOSO'; ...
    'Iterations', noIterations; ...
    'Permutation scheme', 'Cohort-blocked whole-subject profiles'; ...
    'Permutation RNG', sprintf('rng(iteration+%d,''%s'')', ...
        rngSettings.permutation_seed_offset, ...
        rngSettings.loso_permutation_generator); ...
    'Feature selection correlation', ...
        'Pearson for OBN/VRS/MEQ30, Spearman for BDE/GDE/AED'; ...
    'Performance correlation', 'Pearson'; ...
    'Difference behavior fallback', ...
        'prefer *_difference; otherwise raw, which equals the difference'; ...
    'Not assessed', 'OBN was not measured in the LAM cohort'; ...
    'Pooled vectors', ...
        'BDE and OBN have no pooled vector; pooled rows derive from the two single-drug sets'; ...
    'P values', ['Raw in the complete table; Table 2 uses BH across three ' ...
        'outcomes separately within each class']; ...
    'Sweep interpretation', ...
        'Dependent re-analyses of one sample; not independent tests'; ...
    'Workers', workers; ...
    'Expected result rows', expectedRows};
end


function results = load_results_checkpoint(outputMat, header)
if isfile(outputMat)
    loaded = load(outputMat, 'results');
    results = loaded.results;
    assert(iscell(results) && size(results, 2) == numel(header) && ...
        isequal(string(results(1, :)), string(header)), ...
        'CPM:SuppLOSO:CheckpointHeader', ...
        'Existing checkpoint has an unexpected header.');
else
    results = header;
end
end


function tf = analysis_is_complete(results, analysisId, thresholds, config)
if size(results, 1) < 2, tf = false; return; end
ids = string(results(2:end, 1));
testSets = string(results(2:end, 6));
storedThresholds = cell2mat(results(2:end, 13));
expectedTests = [config.drugAStem, config.drugBStem, config.combinedStem];
tf = true;
for threshold = thresholds
    for testSet = expectedTests
        tf = tf && sum(ids == string(analysisId) & ...
            testSets == string(testSet) & ...
            abs(storedThresholds - threshold) < eps) == 1;
    end
end
end


function results = remove_analysis_rows(results, analysisId)
if size(results, 1) < 2, return; end
keep = [true; string(results(2:end, 1)) ~= string(analysisId)];
results = results(keep, :);
end


% ------------------------------------------------------------ checkpointing

function checkpoint_results(results, parameters, outputMat, outputCsv, ...
        outputXlsx, outputTxt)
outputDir = fileparts(outputMat);
temporaryMat = [tempname(outputDir) '.mat'];
save(temporaryMat, 'results', 'parameters', '-v7.3');
atomic_replace(temporaryMat, outputMat, 'MAT checkpoint');

temporaryCsv = [tempname(outputDir) '.csv'];
writecell(results, temporaryCsv);
atomic_replace(temporaryCsv, outputCsv, 'CSV checkpoint');

temporaryXlsx = [tempname(outputDir) '.xlsx'];
writecell(results, temporaryXlsx, 'Sheet', 'Results');
writecell(parameters, temporaryXlsx, 'Sheet', 'Parameters');
atomic_replace(temporaryXlsx, outputXlsx, 'workbook checkpoint');
write_results_text(results, outputTxt);
end


function write_results_text(results, outputTxt)
outputDir = fileparts(outputTxt);
temporaryTxt = [tempname(outputDir) '.txt'];
fileID = fopen(temporaryTxt, 'w');
assert(fileID ~= -1, 'CPM:SuppLOSO:TextOpen', ...
    'Could not open temporary text report.');
cleanupObj = onCleanup(@() fclose_if_open(fileID));
fprintf(fileID, 'Cross-drug LOSO edge-threshold analysis\n');
fprintf(fileID, ['LSD-PCB difference connectome; indicator-coded ' ...
    'covariates; main p<.01 with 0.05 and 0.005 controls.\n']);
fprintf(fileID, ['Final performance uses Pearson r; feature selection ' ...
    'uses the scale-specific correlation.\n']);
for rowIdx = 2:size(results, 1)
    row = results(rowIdx, :);
    fprintf(fileID, '\n------------------------------------------------\n');
    fprintf(fileID, '%s | %s\n', string(row{1}), string(row{2}));
    fprintf(fileID, ['Test=%s; cohort=%s; form=%s; covariates=%s; ' ...
        'selection=%s; edge p<%g; iterations=%d\n'], string(row{6}), ...
        string(row{8}), string(row{9}), string(row{10}), ...
        string(row{11}), row{13}, row{14});
    fprintf(fileID, 'Train n/fold = %d; test n = %d\n', row{15}, row{16});
    fprintf(fileID, 'R_comb = %.4f (p = %.4f)\n', row{18}, row{19});
    fprintf(fileID, 'R_pos  = %.4f (p = %.4f)\n', row{20}, row{21});
    fprintf(fileID, 'R_neg  = %.4f (p = %.4f)\n', row{22}, row{23});
    fprintf(fileID, 'MSE_comb = %.4f (p = %.4f)\n', row{24}, row{25});
    fprintf(fileID, 'MSE_pos  = %.4f (p = %.4f)\n', row{26}, row{27});
    fprintf(fileID, 'MSE_neg  = %.4f (p = %.4f)\n', row{28}, row{29});
    fprintf(fileID, 'MSE_mean-only baseline = %.4f\n', row{30});
end
fclose(fileID);
clear cleanupObj
atomic_replace(temporaryTxt, outputTxt, 'text checkpoint');
end


function atomic_replace(source, target, label)
[ok, message] = movefile(source, target, 'f');
assert(ok, 'CPM:SuppLOSO:CheckpointReplace', ...
    'Could not replace %s: %s', label, message);
end


function save_audit_checkpoint(finalPath, audit)
auditDir = fileparts(finalPath);
temporaryPath = [tempname(auditDir) '.mat'];
save(temporaryPath, 'audit', '-v7.3');
atomic_replace(temporaryPath, finalPath, 'audit checkpoint');
end


% ------------------------------------------------------------------ loading

function [value, variableName] = load_first_var_with_name(pathName)
assert(isfile(pathName), 'CPM:SuppLOSO:MissingInput', ...
    'Missing input: %s', pathName);
loaded = load(pathName);
names = fieldnames(loaded);
assert(isscalar(names), 'CPM:SuppLOSO:AmbiguousInput', ...
    '%s must contain one variable.', pathName);
variableName = names{1};
value = loaded.(variableName);
end


function value = load_first_var(pathName)
[value, ~] = load_first_var_with_name(pathName);
end


function behavior = load_behavior(pathName, behaviorName)
behavior = load_first_var(pathName);
assert(isnumeric(behavior) && isvector(behavior), ...
    'CPM:SuppLOSO:InvalidBehavior', ...
    '%s must be a numeric vector.', behaviorName);
behavior = behavior(:);
assert(all(isfinite(behavior)) && std(behavior) > 0, ...
    'CPM:SuppLOSO:BehaviorValues', ...
    '%s must be finite and nonconstant.', behaviorName);
end


function [behavior, sourceName] = load_behavior_with_fallback( ...
        inputPath, preferredName, fallbackName)
preferredPath = fullfile(inputPath, 'behav', [preferredName '.mat']);
if isfile(preferredPath)
    sourceName = preferredName;
else
    sourceName = fallbackName;
end
behavior = load_behavior(fullfile(inputPath, 'behav', ...
    [sourceName '.mat']), sourceName);
end


function validate_connectome(mats, matrixName, expectedSubjects)
assert(isnumeric(mats) && ndims(mats) == 3 && ...
    size(mats, 1) == size(mats, 2) && ...
    size(mats, 3) == expectedSubjects && all(isfinite(mats(:))), ...
    'CPM:SuppLOSO:InvalidConnectome', ...
    '%s must be a finite square 3-D array with %d subjects.', ...
    matrixName, expectedSubjects);
end


function labels = cohort_labels_from_lsd_rows(globalRows)
labels = strings(numel(globalRows), 1);
labels(globalRows <= 19) = "5HT";
labels(globalRows >= 20 & globalRows <= 42) = "LAM";
labels(globalRows >= 43 & globalRows <= 67) = "LPM";
assert(all(labels ~= ""), 'CPM:SuppLOSO:UnknownLSDRow', ...
    'An LSD row has no cohort label.');
end


function validate_crosswalk_schema(crosswalk)
required = ["test_set", "test_row", "drug", "cohort", ...
    "lsd_row"];
assert(all(ismember(required, ...
    string(crosswalk.Properties.VariableNames))), ...
    'CPM:SuppLOSO:CrosswalkSchema', ...
    'Crosswalk is missing required columns.');
end


function [heldoutRows, anonymousSubjectIndex] = validate_cohort_crosswalk( ...
        crosswalk, config)
rowsA = crosswalk(crosswalk.test_set == config.combinedStem & ...
    lower(crosswalk.drug) == lower(config.drugACrosswalk), :);
rowsB = crosswalk(crosswalk.test_set == config.combinedStem & ...
    lower(crosswalk.drug) == lower(config.drugBCrosswalk), :);
rowsA = sortrows(rowsA, 'test_row');
rowsB = sortrows(rowsB, 'test_row');
n = config.expectedSubjects;
assert(height(rowsA) == n && height(rowsB) == n && ...
    all(rowsA.cohort == config.cohort) && ...
    all(rowsB.cohort == config.cohort), ...
    'CPM:SuppLOSO:CrosswalkCount', ...
    'Crosswalk rows do not match the expected cohort.');
assert(isequal(double(rowsA.test_row(:))', 1:n) && ...
    isequal(double(rowsB.test_row(:))', (n + 1):(2 * n)) && ...
    isequal(rowsA.lsd_row, rowsB.lsd_row) && ...
    numel(unique(rowsA.lsd_row)) == n, ...
    'CPM:SuppLOSO:CrosswalkOrder', ...
    'Crosswalk drug blocks or subject mappings are inconsistent.');
heldoutRows = double(rowsA.lsd_row(:));
anonymousSubjectIndex = (1:n)';
end


function validate_single_drug_connectomes(connectomeDirectory, config, ...
        matrixSuffix, matsA, matsB)
combinedStem = char(config.combinedStem + matrixSuffix);
combinedPath = fullfile(connectomeDirectory, [combinedStem '.mat']);
if ~isfile(combinedPath)
    return % The engine pools the two single-drug arrays directly.
end
combined = load_first_var(combinedPath);
n = config.expectedSubjects;
validate_connectome(combined, combinedStem, 2 * n);
assert(isequaln(matsA, combined(:, :, 1:n)) && ...
    isequaln(matsB, combined(:, :, (n + 1):(2 * n))), ...
    'CPM:SuppLOSO:CombinedConnectomeOrder', ...
    'Single-drug connectomes do not match the combined test-set blocks.');
end


function [combinedSource, crossChecked] = validate_single_drug_behaviors( ...
        inputPath, config, behaviorScale, behaviorA, behaviorB, connectomeForm)
%VALIDATE_SINGLE_DRUG_BEHAVIORS Cross-check the pooled vector when it exists.
% The engine derives the pooled result by concatenating the two single-drug
% test sets, so a stored pooled vector is a consistency check rather than an
% input. BDE and OBN have no pooled vector; the check is then skipped and
% the caller records that on every affected row.
combinedStandard = char(config.combinedStem + "_" + behaviorScale);
if strcmp(connectomeForm,'difference')
    combinedPreferred = [combinedStandard '_difference'];
else
    combinedPreferred = combinedStandard;
end
if isfile(fullfile(inputPath, 'behav', [combinedPreferred '.mat']))
    combinedSource = combinedPreferred;
elseif isfile(fullfile(inputPath, 'behav', [combinedStandard '.mat']))
    combinedSource = combinedStandard;
else
    combinedSource = sprintf('%s (derived from %s and %s)', ...
        combinedStandard, config.drugAStem, config.drugBStem);
    crossChecked = false;
    return
end
combinedBehavior = load_behavior(fullfile(inputPath, 'behav', ...
    [combinedSource '.mat']), combinedSource);
n = config.expectedSubjects;
assert(numel(combinedBehavior) == 2 * n && ...
    isequaln(behaviorA(:), combinedBehavior(1:n)) && ...
    isequaln(behaviorB(:), combinedBehavior((n + 1):(2 * n))), ...
    'CPM:SuppLOSO:CombinedBehaviorOrder', ...
    'Single-drug behaviors do not match the combined test-set blocks.');
crossChecked = true;
end


function fclose_if_open(fileID)
if fileID ~= -1
    try
        fclose(fileID);
    catch
    end
end
end


function value = timestamp_now
value = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS'));
end
