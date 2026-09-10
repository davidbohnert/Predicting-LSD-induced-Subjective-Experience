function consensus_sets = cpm_prepare_consensus_edge_sets( ...
        fold_edge_sets, no_edges, no_iterations, consensus_count)
% CPM_PREPARE_CONSENSUS_EDGE_SETS
% Collapse fold-wise CPM selections into one sliced edge set per iteration.
%
% This performs the same integer fold vote as the submitted external core
% (at least 8 of 10 folds for the paper analyses). Only the storage layout
% changes: full masks are replaced by selected-edge indices.
%
% Keeping the final edge indices by iteration lets a subsequent PARFOR
% receive only the edge set needed for that iteration, rather than
% broadcasting every fold's selections to every worker.

no_folds = numel(fold_edge_sets);
consensus_sets = cell(no_iterations, 1);

for iteration = 1:no_iterations
    pos_votes = zeros(no_edges, 1, 'uint8');
    neg_votes = zeros(no_edges, 1, 'uint8');
    for fold = 1:no_folds
        pos_index = fold_edge_sets{fold}.pos{iteration};
        neg_index = fold_edge_sets{fold}.neg{iteration};
        pos_votes(pos_index) = pos_votes(pos_index) + 1;
        neg_votes(neg_index) = neg_votes(neg_index) + 1;
    end

    consensus_sets{iteration} = struct( ...
        'pos', uint32(find(pos_votes >= consensus_count)), ...
        'neg', uint32(find(neg_votes >= consensus_count)));
end
end
