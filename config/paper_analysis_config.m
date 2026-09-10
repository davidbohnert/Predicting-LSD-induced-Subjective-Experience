function cfg = paper_analysis_config(inputRoot, outputRoot)
%PAPER_ANALYSIS_CONFIG Defaults and presets for the revised manuscript.
%
% CFG = PAPER_ANALYSIS_CONFIG(INPUTROOT, OUTPUTROOT) returns paths, outcome
% metadata, reproducibility settings, and named analysis presets. INPUTROOT
% is the directory containing behav/, connectomes/, covars.mat, and
% covars_indicator.mat. OUTPUTROOT is created by individual pipelines.

configDir = fileparts(mfilename('fullpath'));
repositoryRoot = fileparts(configDir);
if nargin < 1 || isempty(inputRoot)
    inputRoot = fullfile(repositoryRoot, 'data');
end
if nargin < 2 || isempty(outputRoot)
    outputRoot = fullfile(repositoryRoot, 'outputs');
end

cfg.repository_root = repositoryRoot;
cfg.input_root = char(inputRoot);
cfg.output_root = char(outputRoot);
cfg.crosswalk_file = fullfile(repositoryRoot, 'metadata', ...
    'cross_drug_row_map.csv');
cfg.atlas_network_file = fullfile(repositoryRoot, 'metadata', 'atlas', ...
    'Schaefer416_8networks.txt');
cfg.atlas_parcel_metadata_file = fullfile(repositoryRoot, 'metadata', ...
    'atlas', 'atlas_parcel_metadata.csv');
altparcelRoot = fullfile(cfg.input_root, 'connectomes_altparcel');
cfg.altparcel.no_cerebellum_directory = fullfile(altparcelRoot, ...
    'no_cerebellum');
cfg.altparcel.with_cerebellum_directory = fullfile(altparcelRoot, ...
    'with_cerebellum');
cfg.altparcel.with_cerebellum_metadata = fullfile(repositoryRoot, ...
    'metadata', 'atlas', 'Schaefer439_structure_membership.csv');

cfg.runtime.main_iterations = 10000;
cfg.runtime.workers = 5;
cfg.runtime.python_executable = 'python3';
cfg.rng.cv_seed = 123;
cfg.rng.cv_generator = 'twister';
cfg.rng.permutation_seed_offset = 200;
cfg.rng.internal_permutation_generator = 'threefry';
cfg.rng.loso_permutation_generator = 'twister';

cfg.model.folds = 10;
cfg.model.edge_threshold = 0.01;
cfg.model.consensus_fraction = 0.8;
cfg.model.primary_connectome = 'LSD_difference';
cfg.model.covariate_file = fullfile(cfg.input_root, ...
    'covars_indicator.mat');
cfg.model.covariate_columns = ...
    {'study_2','study_3','sex','age','mean_FD'};
cfg.model.grouping_file = fullfile(cfg.input_root, 'covars.mat');

labels = {'BDE','GDE','OBN','AED','VRS','MEQ30'};
corrTypes = {'Spearman','Spearman','Pearson','Spearman','Pearson','Pearson'};
rows = {1:67, 1:67, 1:67, 1:67, 1:67, 20:67};
for index = 1:numel(labels)
    cfg.outcomes(index).label = labels{index};
    cfg.outcomes(index).correlation_type = corrTypes{index};
    cfg.outcomes(index).subject_rows = rows{index};
    cfg.outcomes(index).standard_behavior = ['LSD_' labels{index}];
    cfg.outcomes(index).difference_behavior = ...
        ['LSD_' labels{index} '_difference'];
end

cfg.presets.primary.outcomes = labels;
cfg.presets.primary.connectomes = {cfg.model.primary_connectome};
cfg.presets.primary.adjustment = true;
cfg.presets.primary.thresholds = cfg.model.edge_threshold;
cfg.presets.primary.folds = cfg.model.folds;

cfg.presets.full_internal_grid.outcomes = labels;
cfg.presets.full_internal_grid.connectomes = {'LSD','LSD_difference'};
cfg.presets.full_internal_grid.adjustment = [false true];
cfg.presets.full_internal_grid.thresholds = cfg.model.edge_threshold;
cfg.presets.full_internal_grid.folds = cfg.model.folds;

cfg.presets.fold_sensitivity = cfg.presets.primary;
cfg.presets.fold_sensitivity.folds = [5 10 Inf];
cfg.presets.threshold_sensitivity = cfg.presets.primary;
cfg.presets.threshold_sensitivity.thresholds = [0.005 0.01 0.05];
cfg.presets.alternative_correlation = cfg.presets.primary;
cfg.presets.alternative_correlation.invert_correlation_type = true;
cfg.presets.gsr = cfg.presets.primary;
cfg.presets.gsr.connectomes = {'LSD_GSR_difference'};
cfg.presets.matched_sample = cfg.presets.primary;
cfg.presets.matched_sample.subject_rows = 20:67;
cfg.presets.matched_sample.outcomes = {'GDE','VRS','MEQ30'};
cfg.presets.study_loso = cfg.presets.primary;
cfg.presets.study_loso.fold_source = 'study_labels';
cfg.presets.altparcel_no_cerebellum = cfg.presets.primary;
cfg.presets.altparcel_no_cerebellum.connectomes = ...
    {'LSD_difference_altparcel_no_cerebellum'};
cfg.presets.altparcel_with_cerebellum = cfg.presets.primary;
cfg.presets.altparcel_with_cerebellum.connectomes = ...
    {'LSD_difference_altparcel_with_cerebellum'};

cfg.partition.seeds = 123:222;
cfg.partition.iterations = 1000;
cfg.partition.outcomes = labels;
cfg.learning.repetitions = 100;
cfg.learning.standard_sample_sizes = [20 30 40 50 60 67];
cfg.learning.meq30_sample_sizes = [20 30 40 48];
cfg.learning.iterations = 1000;
cfg.network_interactions.outcomes = {'GDE','MEQ30','VRS'};
cfg.network_interactions.pair_count = 36;
cfg.network_interactions.iterations = 10000;
cfg.network_overlap.outcomes = {'GDE','MEQ30','VRS'};
cfg.network_overlap.thresholds = [0.01 0.05];
cfg.network_overlap.iterations = 10000;
end
