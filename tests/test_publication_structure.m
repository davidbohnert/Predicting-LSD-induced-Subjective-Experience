function test_publication_structure
%TEST_PUBLICATION_STRUCTURE Validate publication defaults and private-path hygiene.
testDir = fileparts(mfilename('fullpath'));
repositoryRoot = fileparts(testDir);
addpath(fullfile(repositoryRoot, 'config'));
cfg = paper_analysis_config('/tmp/cpm_publication_data', ...
    '/tmp/cpm_publication_outputs');
assert(cfg.runtime.main_iterations == 10000 && cfg.runtime.workers == 5);
assert(cfg.rng.cv_seed == 123 && strcmp(cfg.rng.cv_generator, 'twister'));
assert(isequal(cfg.outcomes(end).subject_rows, 20:67));
assert(isequal(cfg.partition.seeds, 123:222));
assert(isequal(cfg.network_overlap.thresholds, [0.01 0.05]));
assert(isfile(cfg.altparcel.with_cerebellum_metadata));
atlasMembership = readtable(cfg.altparcel.with_cerebellum_metadata, ...
    'TextType', 'string');
assert(height(atlasMembership) == 439);
assert(all(ismember({'ParcelIndex','Structure'}, ...
    atlasMembership.Properties.VariableNames)));
assert(numel(unique(atlasMembership.ParcelIndex)) == 439);
assert(nnz(atlasMembership.Structure == "Cerebellum") == 7);
assert(~isfield(cfg.runtime, 'smoke_iterations'));
assert(~isfield(cfg.model, 'primary_network'));

files = [dir(fullfile(repositoryRoot, '**', '*.m')); ...
         dir(fullfile(repositoryRoot, '**', '*.py'))];
for index = 1:numel(files)
    pathName = fullfile(files(index).folder, files(index).name);
    if strcmp(pathName, [mfilename('fullpath') '.m']), continue; end
    text = fileread(pathName);
    assert(~contains(text, '/Users/davidbohnert'), ...
        'CPM:Test:PrivatePath', 'Private path in %s.', pathName);
    assert(~contains(text, 'results/current') && ...
        ~contains(text, 'results/historical'), ...
        'CPM:Test:HistoricalDependency', ...
        'Historical result dependency in %s.', pathName);
end
fprintf('Publication structure and path hygiene passed.\n');
end
