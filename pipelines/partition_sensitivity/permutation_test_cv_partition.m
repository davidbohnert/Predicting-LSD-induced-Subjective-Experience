function [results, details] = permutation_test_cv_partition( ...
        matrix_name, behav_name, all_mats, all_behav, k, thresholds, covars, ...
        output_txt, results, no_iterations, cv_seed, write_text_log, ...
        corr_type, max_workers, rng_settings)
%PERMUTATION_TEST_CV_PARTITION Shared CPM with an explicit partition seed.
if nargin < 10 || isempty(no_iterations), no_iterations = 1000; end
if nargin < 11 || isempty(cv_seed)
    error('CPM:CVSeedRequired', 'Supply the partition seed.');
end
if nargin < 12 || isempty(write_text_log), write_text_log = false; end
if nargin < 13, corr_type = ''; end
if nargin < 14 || isempty(max_workers), max_workers = 0; end
if nargin < 15 || isempty(rng_settings)
    rng_settings = struct('cv_generator','twister', ...
        'permutation_seed_offset',200,'internal_permutation_generator','threefry');
end
validateattributes(cv_seed, {'numeric'}, {'scalar','integer','nonnegative','finite'});
rng_settings.cv_seed = cv_seed;
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root, 'scripts', 'common'));
[results, details] = cpm_permutation_cv(matrix_name, behav_name, ...
    all_mats, all_behav, k, thresholds, covars, output_txt, results, ...
    no_iterations, corr_type, max_workers, rng_settings, write_text_log);
end
