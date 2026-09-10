function [results, details] = cpm_permutation_cv( ...
    matrix_name, behav_name, all_mats, all_behav, k, thresholds, covars, ...
    output_txt, results, no_iterations, corr_type_override, max_workers, ...
    rng_settings, write_text_log)
%CPM_PERMUTATION_CV Fixed-partition CPM with cached edge correlations.
% Select edges, summarize network strengths, fit the three linear models,
% and predict held-out subjects. Each threshold has separate masks and fits.
% Column 1 is observed; remaining columns are behavior permutations.

if nargin < 14, write_text_log = true; end
if nargin < 10 || isempty(no_iterations)
    no_iterations = 1000;
end
if nargin < 11
    corr_type_override = '';
end
if nargin < 12 || isempty(max_workers), max_workers = 0; end
if nargin < 13 || isempty(rng_settings)
    rng_settings = default_rng_settings();
end
validateattributes(max_workers, {'numeric'}, ...
    {'scalar', 'integer', '>=', 0});
validateattributes(no_iterations, {'numeric'}, ...
    {'scalar', 'integer', '>=', 2}, mfilename, 'no_iterations');
validateattributes(thresholds, {'numeric'}, ...
    {'vector', 'real', 'finite', '>', 0, '<', 1}, ...
    mfilename, 'thresholds');
thresholds = thresholds(:)';

all_behav = all_behav(:);
no_sub = size(all_mats, 3);
if numel(all_behav) ~= no_sub || ...
        (~isempty(covars) && size(covars, 1) ~= no_sub)
    error('CPM:InputSizeMismatch', ...
        'Connectivity, behavior, and covariates must contain the same subjects.');
end
if any(~isfinite(all_mats), 'all') || any(~isfinite(all_behav), 'all') || ...
        (~isempty(covars) && any(~isfinite(covars), 'all'))
    error('CPM:NonfiniteInput', 'CPM inputs must be finite.');
end

repository_dir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repository_dir, 'scripts', 'common'));

%% 1. FIXED PARTITION, CORRELATION TYPE, AND PERMUTATIONS
rng(rng_settings.cv_seed, rng_settings.cv_generator);
cv_fixed = cvpartition(no_sub, 'KFold', k);

covar_text = '';
if ~isempty(covars), covar_text = ' - with covariates'; end
if ~isempty(corr_type_override)
    corr_type = char(validatestring(corr_type_override, ...
        {'Pearson', 'Spearman'}));
elseif contains(behav_name, 'LSD_OBN') || ...
        contains(behav_name, 'LSD_MEQ30') || ...
        contains(behav_name, 'LSD_VRS')
    corr_type = 'Pearson';
else
    corr_type = 'Spearman';
end

fprintf(['Running Permutation: %s vs %s%s ' ...
    '(k=%d, thresholds=%s, iterations=%d)\n'], ...
    matrix_name, behav_name, covar_text, k, ...
    mat2str(thresholds), no_iterations);

perm_behav = cpm_behavior_permutations( ...
    all_behav, no_iterations, rng_settings);
edge_layout = cpm_prepare_edge_layout(all_mats, corr_type);
all_correlation_edges = edge_layout.correlation_edges;
all_strength_edges = edge_layout.strength_edges;
strength_divisor = edge_layout.strength_divisor;
no_thresholds = numel(thresholds);

%% 2. SHEN-STYLE CPM WITH ONE FIXED-SIDE CACHE PER FOLD
fold_output = cell(k, 1);
parfor (fold = 1:k, max_workers)
    train_idx = find(training(cv_fixed, fold));
    test_idx = find(test(cv_fixed, fold));
    train_correlation_edges = ...
        all_correlation_edges(:, train_idx)';
    train_strength_edges = all_strength_edges(:, train_idx);
    test_strength_edges = all_strength_edges(:, test_idx);

    if isempty(covars)
        train_covars = [];
    else
        train_covars = covars(train_idx, :);
    end
    correlation_cache = cpm_prepare_correlation_cache( ...
        train_correlation_edges, train_covars, corr_type);

    fold_output{fold} = run_one_fold( ...
        train_strength_edges, test_strength_edges, ...
        strength_divisor, train_idx, test_idx, ...
        perm_behav, thresholds, correlation_cache);
    fprintf('  ...Cached fold %d / %d complete\n', fold, k);
end

%% 3. REASSEMBLE PREDICTIONS AND CALCULATE PERFORMANCE
behav_pred_comb = zeros(no_sub, no_iterations, no_thresholds);
behav_pred_pos = zeros(no_sub, no_iterations, no_thresholds);
behav_pred_neg = zeros(no_sub, no_iterations, no_thresholds);
for fold = 1:k
    test_idx = fold_output{fold}.test_idx;
    behav_pred_comb(test_idx, :, :) = fold_output{fold}.comb;
    behav_pred_pos(test_idx, :, :) = fold_output{fold}.pos;
    behav_pred_neg(test_idx, :, :) = fold_output{fold}.neg;
end

perm_r = zeros(no_iterations, 3, no_thresholds);
perm_mse = zeros(no_iterations, 3, no_thresholds);
for threshold_index = 1:no_thresholds
    for iteration = 1:no_iterations
        current_behavior = perm_behav(:, iteration);
        pred_comb = behav_pred_comb(:, iteration, threshold_index);
        pred_pos = behav_pred_pos(:, iteration, threshold_index);
        pred_neg = behav_pred_neg(:, iteration, threshold_index);
        perm_r(iteration, :, threshold_index) = [ ...
            corr(pred_comb, current_behavior), ...
            corr(pred_pos, current_behavior), ...
            corr(pred_neg, current_behavior)];
        perm_mse(iteration, :, threshold_index) = [ ...
            mean((current_behavior - pred_comb).^2), ...
            mean((current_behavior - pred_pos).^2), ...
            mean((current_behavior - pred_neg).^2)];
    end
end

%% 4. APPEND ONE CONVENTIONAL RESULT ROW PER THRESHOLD
for threshold_index = 1:no_thresholds
    current_r = perm_r(:, :, threshold_index);
    current_mse = perm_mse(:, :, threshold_index);
    true_r = current_r(1, :);
    true_mse = current_mse(1, :);
    p_r = sum(current_r >= true_r, 1) / no_iterations;
    p_mse = sum(current_mse <= true_mse, 1) / no_iterations;
    threshold = thresholds(threshold_index);

    if write_text_log && ~isempty(output_txt)
        append_text_result(output_txt, matrix_name, behav_name, covar_text, ...
            k, threshold, true_r, p_r, true_mse, p_mse);
    end
    results(end+1, :) = {sprintf('%s vs %s%s', ...
        matrix_name, behav_name, covar_text), k, threshold, ...
        true_r(1), p_r(1), true_r(2), p_r(2), true_r(3), p_r(3), ...
        true_mse(1), p_mse(1), true_mse(2), p_mse(2), ...
        true_mse(3), p_mse(3)}; %#ok<AGROW>
end

if nargout > 1
    details = struct( ...
        'perm_r', squeeze_scalar_threshold(perm_r, no_thresholds), ...
        'perm_mse', squeeze_scalar_threshold(perm_mse, no_thresholds), ...
        'perm_behav', perm_behav, ...
        'observed_behavior', all_behav, ...
        'observed_predictions', struct( ...
            'combined', reshape(behav_pred_comb(:,1,:),no_sub,no_thresholds), ...
            'positive', reshape(behav_pred_pos(:,1,:),no_sub,no_thresholds), ...
            'negative', reshape(behav_pred_neg(:,1,:),no_sub,no_thresholds)), ...
        'corr_type', corr_type, ...
        'cv_partition', cv_fixed, ...
        'cv_seed', rng_settings.cv_seed, ...
        'thresholds', thresholds, ...
        'uses_unique_edges', edge_layout.uses_unique_edges);
end
end


function output = run_one_fold(train_edges, test_edges, strength_divisor, ...
        train_idx, test_idx, perm_behav, thresholds, correlation_cache)
num_test = numel(test_idx);
no_iterations = size(perm_behav, 2);
no_thresholds = numel(thresholds);
saved_warning = warning('off', 'stats:regress:RankDefDesignMat');
warning_cleanup = onCleanup(@() warning(saved_warning));

pred_comb = zeros(num_test, no_iterations, no_thresholds);
pred_pos = zeros(num_test, no_iterations, no_thresholds);
pred_neg = zeros(num_test, no_iterations, no_thresholds);

for iteration = 1:no_iterations
    train_behavior = perm_behav(train_idx, iteration);

    % Step 1: select positive and negative edges.
    [~, pos_masks, neg_masks] = cpm_cached_edge_selection( ...
        correlation_cache, train_behavior, thresholds);

    for threshold_index = 1:no_thresholds
        pos_mask = pos_masks{threshold_index};
        neg_mask = neg_masks{threshold_index};

        % Step 2: summarize each network from only its selected edges.
        train_sumpos = ...
            sum(train_edges(pos_mask, :), 1)' / strength_divisor;
        train_sumneg = ...
            sum(train_edges(neg_mask, :), 1)' / strength_divisor;
        test_sumpos = ...
            sum(test_edges(pos_mask, :), 1)' / strength_divisor;
        test_sumneg = ...
            sum(test_edges(neg_mask, :), 1)' / strength_divisor;

        % Step 3: fit the conventional combined, positive, and negative models.
        [b_comb, b_pos, b_neg] = cpm_fit_three_models( ...
            train_behavior, train_sumpos, train_sumneg);

        % Step 4: predict held-out subjects.
        pred_comb(:, iteration, threshold_index) = ...
            b_comb(1) * test_sumpos + b_comb(2) * test_sumneg + b_comb(3);
        pred_pos(:, iteration, threshold_index) = ...
            b_pos(1) * test_sumpos + b_pos(2);
        pred_neg(:, iteration, threshold_index) = ...
            b_neg(1) * test_sumneg + b_neg(2);
    end
end

output = struct('test_idx', test_idx, ...
    'comb', pred_comb, 'pos', pred_pos, 'neg', pred_neg);
clear warning_cleanup
end


function append_text_result(output_txt, matrix_name, behav_name, ...
        covar_text, k, threshold, true_r, p_r, true_mse, p_mse)
try
    fileID = fopen(output_txt, 'a');
    if fileID < 0, error('Could not open output file: %s', output_txt); end
    cleanup = onCleanup(@() fclose(fileID));
    fprintf(fileID, '\n------------------------------------------------\n');
    fprintf(fileID, '%s vs %s%s\nParams: k=%d, THR=%.3f\n', ...
        matrix_name, behav_name, covar_text, k, threshold);
    fprintf(fileID, 'R_comb = %1.4f (p = %1.4f)\n', true_r(1), p_r(1));
    fprintf(fileID, 'R_pos  = %1.4f (p = %1.4f)\n', true_r(2), p_r(2));
    fprintf(fileID, 'R_neg  = %1.4f (p = %1.4f)\n', true_r(3), p_r(3));
    fprintf(fileID, 'MSE_comb = %1.4f (p = %1.4f)\n', ...
        true_mse(1), p_mse(1));
    fprintf(fileID, 'MSE_pos  = %1.4f (p = %1.4f)\n', ...
        true_mse(2), p_mse(2));
    fprintf(fileID, 'MSE_neg  = %1.4f (p = %1.4f)\n', ...
        true_mse(3), p_mse(3));
    clear cleanup
catch ME
    warning('CPM:LogWriteFailed', ...
        'Failed to write to text file: %s', ME.message);
end
end


function value = squeeze_scalar_threshold(value, no_thresholds)
if no_thresholds == 1
    value = value(:, :, 1);
end
end


function settings = default_rng_settings()
settings = struct('cv_seed',123, 'cv_generator','twister', ...
    'permutation_seed_offset',200, ...
    'internal_permutation_generator','threefry');
end
