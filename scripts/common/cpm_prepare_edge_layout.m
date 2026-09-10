function layout = cpm_prepare_edge_layout(mats, corr_type)
% CPM_PREPARE_EDGE_LAYOUT
% Vectorize connectomes for cached feature selection and network strengths.
%
% The submitted implementation correlated every entry of each symmetric
% node-by-node matrix, built a full symmetric mask, summed both triangles,
% and divided by two. For validated symmetric inputs, one off-diagonal
% triangle contains the same undirected edges and its single sum is the same
% network-strength definition.
%
% Numerically symmetric Pearson matrices and Spearman matrices with
% identical mirrored rank sequences use one unique off-diagonal triangle.
% Full symmetric masks are reconstructed for Shen-style output. Network
% strengths use the mean of each mirrored pair, representing the original
% full-mask sum divided by two apart from floating-point summation order.
% A full-edge fallback protects unusual asymmetric inputs.

corr_type = char(validatestring(corr_type, {'Pearson', 'Spearman'}));
no_node = size(mats, 1);
no_subjects = size(mats, 3);
full_edges = reshape(mats, [], no_subjects);

layout.no_node = no_node;
layout.corr_type = corr_type;
[row_index, column_index] = find(triu(true(no_node), 1));
layout.upper_index = sub2ind( ...
    [no_node, no_node], row_index, column_index);
layout.lower_index = sub2ind( ...
    [no_node, no_node], column_index, row_index);
upper_edges = full_edges(layout.upper_index, :);
lower_edges = full_edges(layout.lower_index, :);
diagonal_edges = full_edges(1:(no_node + 1):(no_node * no_node), :);

diagonal_is_constant = all( ...
    diagonal_edges == diagonal_edges(:, 1), 'all');
% Constant diagonal entries could never pass correlation-based selection in
% the submitted full-matrix calculation, so they can be omitted. If a new
% input violates that assumption, retain the full representation.
if strcmpi(corr_type, 'Spearman')
    mirrors_are_equivalent = mirrored_ranks_match( ...
        upper_edges, lower_edges);
else
    mirrors_are_equivalent = numerically_symmetric( ...
        upper_edges, lower_edges);
end
layout.uses_unique_edges = ...
    diagonal_is_constant && mirrors_are_equivalent;

if layout.uses_unique_edges
    symmetric_edges = (upper_edges + lower_edges) / 2;
    if strcmpi(corr_type, 'Spearman')
        % Spearman depends on order rather than numerical closeness. Either
        % mirror has been verified to have exactly the same ranks; retain
        % one side to avoid averaging extremely close values into a new tie.
        layout.correlation_edges = upper_edges;
    else
        layout.correlation_edges = symmetric_edges;
    end
    layout.strength_edges = symmetric_edges;
    layout.strength_divisor = 1;
else
    warning('CPM:FullEdgeFallback', ...
        ['Connectomes are not symmetric under the %s edge criterion; ' ...
         'retaining the full matrix representation.'], corr_type);
    layout.correlation_edges = full_edges;
    layout.strength_edges = full_edges;
    layout.strength_divisor = 2;
end

layout.no_edges = size(layout.correlation_edges, 1);
end


function matches = mirrored_ranks_match(upper_edges, lower_edges)
% Compare ranks in blocks to limit temporary memory.
no_edges = size(upper_edges, 1);
block_size = 4096;
matches = true;
for first_edge = 1:block_size:no_edges
    last_edge = min(first_edge + block_size - 1, no_edges);
    edge_block = first_edge:last_edge;
    upper_ranks = tiedrank(upper_edges(edge_block, :)');
    lower_ranks = tiedrank(lower_edges(edge_block, :)');
    if ~isequal(upper_ranks, lower_ranks)
        matches = false;
        return
    end
end
end


function symmetric = numerically_symmetric(upper_edges, lower_edges)
% Pearson changes smoothly under machine-precision mirror differences. The
% tolerance accepts numerical roundoff but rejects substantively asymmetric
% matrices, which take the full-edge fallback above.
data_scale = max([1, max(abs(upper_edges), [], 'all'), ...
    max(abs(lower_edges), [], 'all')]);
tolerance = 64 * eps(class(upper_edges)) * data_scale;
no_edges = size(upper_edges, 1);
block_size = 4096;
symmetric = true;
for first_edge = 1:block_size:no_edges
    last_edge = min(first_edge + block_size - 1, no_edges);
    edge_block = first_edge:last_edge;
    if any(abs(upper_edges(edge_block, :) - ...
            lower_edges(edge_block, :)) > tolerance, 'all')
        symmetric = false;
        return
    end
end
end
