function test_checkpoint_workflow
%TEST_CHECKPOINT_WORKFLOW Small synthetic resume and input-contract checks.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(root));
scratch = tempname;
mkdir(scratch);
cleanup = onCleanup(@() rmdir(scratch,'s'));
data = fullfile(scratch,'data');
mkdir(data); mkdir(fullfile(data,'behav')); mkdir(fullfile(data,'connectomes'));
rng(71,'twister');
studies = [ones(19,1);2*ones(23,1);3*ones(25,1)];
covars = [double(studies==2),double(studies==3), ...
    double(mod((1:67)',2)),20+rand(67,1)*15,rand(67,1)*.1];
save(fullfile(data,'covars_indicator.mat'),'covars');
covars = [studies,covars(:,3:5)];
save(fullfile(data,'covars.mat'),'covars');
stems = {'LSD','amphetamine','MDMA','psilocybin','mescaline'};
sizes = [67 23 23 25 25];
for i = 1:numel(stems)
    mats = zeros(8,8,sizes(i));
    for j = 1:sizes(i)
        a = randn(8); a = (a+a')/2; a(1:9:end) = 0;
        mats(:,:,j) = a;
    end
    save(fullfile(data,'connectomes',[stems{i} '_difference.mat']),'mats');
    behavior = randn(sizes(i),1);
    save(fullfile(data,'behav',[stems{i} '_GDE_difference.mat']),'behavior');
end
behavior = randn(67,1);
save(fullfile(data,'behav','LSD_VRS_difference.mat'),'behavior');
for label = {'BDE','AED','OBN','MEQ30'}
    if strcmp(label{1},'MEQ30'), n = 48; else, n = 67; end
    behavior = randn(n,1);
    save(fullfile(data,'behav',['LSD_' label{1} '_difference.mat']),'behavior');
end
cfg = paper_analysis_config(data,fullfile(scratch,'outputs'));
cfg.runtime.workers = 0;
cfg.presets.matched_sample.outcomes = {'GDE','VRS'};
cfg.presets.matched_sample.thresholds = .2;
cfg.presets.matched_sample.folds = 3;
out = fullfile(cfg.output_root,'internal_validation','matched_sample');
first = internal_validation_main(cfg,"matched_sample", ...
    no_iterations=2,workers=0,use_parallel=false);
assert(height(first)==2);
% Simulate an interruption after the first completed analysis unit.
path = fullfile(out,'internal_validation_checkpoint.mat');
saved = load(path);
saved.resultRows = saved.resultRows(1,:);
saved.details = saved.details(1);
saved.completedUnits = saved.completedUnits(1);
cpm_atomic_save(path,saved);
second = internal_validation_main(cfg,"matched_sample", ...
    no_iterations=2,workers=5,use_parallel=false);
assert(isequaln(first,second), 'Resume changed or duplicated results.');
expect_error(@() internal_validation_main(cfg,"matched_sample", ...
    no_iterations=3,workers=0,use_parallel=false), 'CPM:Checkpoint:Incompatible');
changed = cfg; changed.rng.permutation_seed_offset = 201;
expect_error(@() internal_validation_main(changed,"matched_sample", ...
    no_iterations=2,workers=0,use_parallel=false), 'CPM:Checkpoint:Incompatible');
behavior = randn(67,1);
save(fullfile(data,'behav','LSD_VRS_difference.mat'),'behavior');
expect_error(@() internal_validation_main(cfg,"matched_sample", ...
    no_iterations=2,workers=0,use_parallel=false), 'CPM:Checkpoint:Incompatible');
% No pooled connectomes/behaviors are present.
cross = cross_drug_loso_main(cfg,no_iterations=2,workers=0, ...
    scales={'GDE'},thresholds=.2);
assert(size(cross,1)==7, 'Expected two classes x three test sets.');
again = cross_drug_loso_main(cfg,no_iterations=2,workers=0, ...
    scales={'GDE'},thresholds=.2);
assert(isequaln(cross,again));
expect_error(@() cross_drug_loso_main(cfg,no_iterations=3,workers=0, ...
    scales={'GDE'},thresholds=.2), 'CPM:Checkpoint:Incompatible');
cfg.partition.outcomes = {'GDE','VRS'};
partition = partition_sensitivity_main(cfg,no_iterations=2,workers=0,seeds=123);
assert(height(partition)==2);
expect_error(@() partition_sensitivity_main(cfg,no_iterations=3,workers=0,seeds=123), ...
    'CPM:Checkpoint:Incompatible');
% The learning-curve entry point requires the paper's 416-node input shape.
% Repeat the small synthetic blocks; use only one 20-person sample and two
% iterations per outcome. This is an interface check, not an analysis rerun.
matrixPath = fullfile(data,'connectomes','LSD_difference.mat');
small = load(matrixPath,'mats');
mats = repmat(small.mats,52,52,1);
save(matrixPath,'mats');
curve = learning_curve_main(cfg,no_iterations=2,workers=0,repetitions=1, ...
    standard_sample_sizes=20,meq30_sample_sizes=20);
assert(height(curve)==6);
expect_error(@() learning_curve_main(cfg,no_iterations=3,workers=0,repetitions=1, ...
    standard_sample_sizes=20,meq30_sample_sizes=20), 'CPM:Checkpoint:Incompatible');
% A prior output folder lacking a signature must not be adopted.
legacy = fullfile(scratch,'legacy'); mkdir(legacy);
save(fullfile(legacy,'old.mat'),'behavior');
expect_error(@() cpm_checkpoint_guard(legacy,struct,{},{}), ...
    'CPM:Checkpoint:Unverified');
clear cleanup
fprintf('Checkpoint resume, incompatible inputs/settings, and unpooled LOSO passed.\n');
end

function expect_error(action, identifier)
try
    action();
catch problem
    assert(strcmp(problem.identifier,identifier), '%s: %s', ...
        problem.identifier,problem.message);
    return
end
error('Expected error %s.',identifier);
end
