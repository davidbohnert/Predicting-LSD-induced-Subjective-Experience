function consensus = cpm_build_consensus_masks( ...
        mats_train, behavior_columns, covars, corr_type, ...
        fixed_partition, p_threshold, consensus_fraction, max_workers)
% CPM_BUILD_CONSENSUS_MASKS
% Build observed or permuted CPM consensus masks from training data only.
%
% BEHAVIOR_COLUMNS is subjects-by-models. A standalone visualization export
% passes one observed behavior vector; permutation analyses pass the
% observed vector followed by shuffled vectors. This keeps both workflows
% on exactly the same feature-selection implementation.
%
% MAX_WORKERS controls only the fold loop:
%   0     - run serially without opening a parallel pool
%   []    - use MATLAB's ordinary PARFOR behavior
%   n > 0 - use at most n workers

if nargin < 8
    max_workers = [];
end

corr_type = char(validatestring(corr_type, {'Pearson', 'Spearman'}));
validateattributes(mats_train, {'numeric'}, ...
    {'3d', 'real', 'finite', 'nonempty'});
validateattributes(behavior_columns, {'numeric'}, ...
    {'2d', 'real', 'finite', 'nonempty'});
validateattributes(p_threshold, {'numeric'}, ...
    {'scalar', 'real', 'finite', '>', 0, '<', 1});
validateattributes(consensus_fraction, {'numeric'}, ...
    {'scalar', 'real', 'finite', '>', 0, '<=', 1});
if ~isempty(max_workers)
    validateattributes(max_workers, {'numeric'}, ...
        {'scalar', 'integer', '>=', 0});
end

no_subjects = size(mats_train, 3);
if isvector(behavior_columns)
    behavior_columns = behavior_columns(:);
end
assert(size(behavior_columns, 1) == no_subjects, ...
    'CPM:InputSizeMismatch', ...
    'Training connectomes and behavior columns do not agree.');
assert(isempty(covars) || ...
    (isnumeric(covars) && ismatrix(covars) && ...
    size(covars, 1) == no_subjects && all(isfinite(covars(:)))), ...
    'CPM:InputSizeMismatch', ...
    'Covariates must be finite and have one row per training subject.');
assert(fixed_partition.NumObservations == no_subjects, ...
    'CPM:PartitionSizeMismatch', ...
    'The fixed partition does not match the training subject count.');

layout = cpm_prepare_edge_layout(mats_train, corr_type);
correlation_edges = layout.correlation_edges;
no_models = size(behavior_columns, 2);
no_folds = fixed_partition.NumTestSets;
fold_edge_sets = cell(no_folds, 1);

if isequal(max_workers, 0)
    for fold = 1:no_folds
        fold_edge_sets{fold} = select_fold_edges( ...
            fold, fixed_partition, correlation_edges, behavior_columns, ...
            covars, corr_type, p_threshold);
        fprintf('  ...Consensus fold %d / %d complete\n', fold, no_folds);
    end
elseif isempty(max_workers)
    parfor fold = 1:no_folds
        fold_edge_sets{fold} = select_fold_edges( ...
            fold, fixed_partition, correlation_edges, behavior_columns, ...
            covars, corr_type, p_threshold);
        fprintf('  ...Consensus fold %d / %d complete\n', fold, no_folds);
    end
else
    parfor (fold = 1:no_folds, max_workers)
        fold_edge_sets{fold} = select_fold_edges( ...
            fold, fixed_partition, correlation_edges, behavior_columns, ...
            covars, corr_type, p_threshold);
        fprintf('  ...Consensus fold %d / %d complete\n', fold, no_folds);
    end
end

consensus_count = ceil(consensus_fraction * no_folds);
edge_sets = cpm_prepare_consensus_edge_sets( ...
    fold_edge_sets, layout.no_edges, no_models, consensus_count);
[observed_pos, observed_neg] = masks_from_indices( ...
    edge_sets{1}, layout.no_edges);

consensus.edge_sets = edge_sets;
consensus.layout = layout;
consensus.observed_pos_mask = ...
    cpm_expand_edge_mask(observed_pos, layout);
consensus.observed_neg_mask = ...
    cpm_expand_edge_mask(observed_neg, layout);
consensus.no_folds = no_folds;
consensus.consensus_count = consensus_count;
consensus.p_threshold = p_threshold;
consensus.consensus_fraction = consensus_fraction;
consensus.corr_type = corr_type;
end


function fold_result = select_fold_edges( ...
        fold, fixed_partition, correlation_edges, behavior_columns, ...
        covars, corr_type, p_threshold)
train_idx = find(training(fixed_partition, fold));
fold_edges = correlation_edges(:, train_idx)';
if isempty(covars)
    fold_covars = [];
else
    fold_covars = covars(train_idx, :);
end
cache = cpm_prepare_correlation_cache( ...
    fold_edges, fold_covars, corr_type);

no_models = size(behavior_columns, 2);
pos_sets = cell(1, no_models);
neg_sets = cell(1, no_models);
for model_index = 1:no_models
    [~, pos_masks, neg_masks] = cpm_cached_edge_selection( ...
        cache, behavior_columns(train_idx, model_index), p_threshold);
    pos_sets{model_index} = uint32(find(pos_masks{1}));
    neg_sets{model_index} = uint32(find(neg_masks{1}));
end
fold_result = struct('pos', {pos_sets}, 'neg', {neg_sets});
end


function [pos_mask, neg_mask] = masks_from_indices(edge_set, no_edges)
pos_mask = false(no_edges, 1);
neg_mask = false(no_edges, 1);
pos_mask(edge_set.pos) = true;
neg_mask(edge_set.neg) = true;
end
