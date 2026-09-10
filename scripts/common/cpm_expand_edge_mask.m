function full_mask = cpm_expand_edge_mask(edge_mask, layout)
% CPM_EXPAND_EDGE_MASK
% Convert a cached edge vector to a full Shen-style symmetric node mask.
%
% Exported masks therefore retain the dimensions and interpretation of the
% full masks written by the submitted external-validation pipeline.

if layout.uses_unique_edges
    full_mask = false(layout.no_node, layout.no_node);
    selected_upper = layout.upper_index(edge_mask);
    selected_lower = layout.lower_index(edge_mask);
    full_mask(selected_upper) = true;
    full_mask(selected_lower) = true;
else
    full_mask = reshape(edge_mask, layout.no_node, layout.no_node);
end
end
