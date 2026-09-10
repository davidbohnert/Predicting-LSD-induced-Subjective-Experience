function result = cpm_build_native_consensus_sets( ...
        matsTrain, behaviorColumns, covars, corrType, fixedPartition, ...
        thresholds, consensusFraction, maxWorkers)
%CPM_BUILD_NATIVE_CONSENSUS_SETS Build ordinary CPM consensus masks.
%
% Unlike the superseded size-matched overlap helper, this function returns
% only edges selected in at least the requested number of folds. It never
% ranks, pads, fills, truncates, or discards a mask.

if nargin < 8
    maxWorkers = [];
end

corrType = char(validatestring(corrType, {'Pearson', 'Spearman'}));
thresholds = thresholds(:)';
layout = cpm_prepare_edge_layout(matsTrain, corrType);
correlationEdges = layout.correlation_edges;
noModels = size(behaviorColumns, 2);
noFolds = fixedPartition.NumTestSets;
noThresholds = numel(thresholds);
foldEdgeSets = cell(noFolds, 1);

if isequal(maxWorkers, 0)
    for fold = 1:noFolds
        foldEdgeSets{fold} = select_fold_edges_multi( ...
            fold, fixedPartition, correlationEdges, behaviorColumns, ...
            covars, corrType, thresholds);
        fprintf('  ...Native-mask fold %d / %d complete\n', fold, noFolds);
    end
elseif isempty(maxWorkers)
    parfor fold = 1:noFolds
        foldEdgeSets{fold} = select_fold_edges_multi( ...
            fold, fixedPartition, correlationEdges, behaviorColumns, ...
            covars, corrType, thresholds);
        fprintf('  ...Native-mask fold %d / %d complete\n', fold, noFolds);
    end
else
    parfor (fold = 1:noFolds, maxWorkers)
        foldEdgeSets{fold} = select_fold_edges_multi( ...
            fold, fixedPartition, correlationEdges, behaviorColumns, ...
            covars, corrType, thresholds);
        fprintf('  ...Native-mask fold %d / %d complete\n', fold, noFolds);
    end
end

consensusCount = ceil(consensusFraction * noFolds);
positiveSets = cell(noModels, noThresholds);
positiveCounts = zeros(noModels, noThresholds, 'uint32');
for modelIndex = 1:noModels
    for thresholdIndex = 1:noThresholds
        votes = zeros(layout.no_edges, 1, 'uint8');
        for fold = 1:noFolds
            selected = foldEdgeSets{fold}.pos{modelIndex, thresholdIndex};
            votes(selected) = votes(selected) + 1;
        end
        indices = uint32(find(votes >= consensusCount));
        positiveSets{modelIndex, thresholdIndex} = indices;
        positiveCounts(modelIndex, thresholdIndex) = uint32(numel(indices));
    end
end

result = struct;
result.layout = layout;
result.thresholds = thresholds;
result.consensus_count = consensusCount;
result.positive_sets = positiveSets;
result.positive_counts = positiveCounts;
result.no_models = noModels;
result.no_folds = noFolds;
result.mask_policy = 'native >=8/10 consensus; no ranking, padding, or filling';
end


function foldResult = select_fold_edges_multi( ...
        fold, fixedPartition, correlationEdges, behaviorColumns, ...
        covars, corrType, thresholds)
trainIdx = find(training(fixedPartition, fold));
foldEdges = correlationEdges(:, trainIdx)';
if isempty(covars)
    foldCovars = [];
else
    foldCovars = covars(trainIdx, :);
end
cache = cpm_prepare_correlation_cache(foldEdges, foldCovars, corrType);

noModels = size(behaviorColumns, 2);
noThresholds = numel(thresholds);
posSets = cell(noModels, noThresholds);
for modelIndex = 1:noModels
    [~, posMasks] = cpm_cached_edge_selection( ...
        cache, behaviorColumns(trainIdx, modelIndex), thresholds);
    for thresholdIndex = 1:noThresholds
        posSets{modelIndex, thresholdIndex} = ...
            uint32(find(posMasks{thresholdIndex}));
    end
end
foldResult = struct('pos', {posSets});
end
