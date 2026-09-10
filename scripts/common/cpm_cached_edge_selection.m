function [edge_r, pos_by_threshold, neg_by_threshold] = ...
        cpm_cached_edge_selection(cache, behavior, thresholds)
% CPM_CACHED_EDGE_SELECTION
% Select CPM edges from a fold-specific fixed-side correlation cache.
%
% This returns the same p<THRESHOLD decision without constructing a full
% edge-p-value array for every permutation. The rare completely untied
% Spearman case retains MATLAB's AS89 approximation.
%
% The submitted code formed p for every edge and then tested p<THRESHOLD.
% For Pearson-style tests that two-sided p cutoff is exactly equivalent to
% abs(r)>r_critical at fixed degrees of freedom. The strict inequality is
% retained so an edge exactly on the boundary is treated the same way.

behavior = behavior(:);
thresholds = thresholds(:)';
validateattributes(thresholds, {'numeric'}, ...
    {'vector', 'real', 'finite', '>', 0, '<', 1});

no_thresholds = numel(thresholds);
pos_by_threshold = cell(1, no_thresholds);
neg_by_threshold = cell(1, no_thresholds);

if strcmpi(cache.corr_type, 'Spearman') && cache.n < 10
    if isempty(cache.covars)
        [edge_r, edge_p] = corr( ...
            cache.edge_data, behavior, 'type', 'Spearman');
    else
        [edge_r, edge_p] = partialcorr( ...
            cache.edge_data, behavior, cache.covars, 'type', 'Spearman');
    end
    for threshold_index = 1:no_thresholds
        selected = edge_p < thresholds(threshold_index);
        pos_by_threshold{threshold_index} = edge_r > 0 & selected;
        neg_by_threshold{threshold_index} = edge_r < 0 & selected;
    end
    return
end

is_partial = ~isempty(cache.covars);
if is_partial || strcmpi(cache.corr_type, 'Pearson')
    edge_r = cached_linear_correlation(cache, behavior);
    df = max(cache.n - cache.covar_rank - 2, 0);
    for threshold_index = 1:no_thresholds
        r_critical = correlation_r_critical( ...
            thresholds(threshold_index), df);
        selected = abs(edge_r) > r_critical;
        pos_by_threshold{threshold_index} = edge_r > 0 & selected;
        neg_by_threshold{threshold_index} = edge_r < 0 & selected;
    end
    return
end

% Ordinary Spearman: calculate the exact rank-distance statistic once.
[behavior_rank, behavior_tieadj] = tiedrank(behavior);
n = cache.n;
n3const = (n + 1) * n * (n - 1) / 3;
D = cache.edge_rank_sumsq + sum(behavior_rank.^2) - ...
    2 * (cache.edge_ranks' * behavior_rank);
meanD = (n3const - (cache.edge_tieadj + behavior_tieadj) / 3) / 2;
stdD = sqrt((n3const / 2 - cache.edge_tieadj / 3) .* ...
    (n3const / 2 - behavior_tieadj / 3) / (n - 1));

constant_adj = (n + 1) * n * (n - 1) / 2;
stdD((cache.edge_tieadj == constant_adj) | ...
    (behavior_tieadj == constant_adj)) = 0;
edge_r = (meanD - D) ./ (sqrt(n - 1) .* stdD);
outside = abs(edge_r) > 1;
edge_r(outside) = sign(edge_r(outside));

if behavior_tieadj > 0
    % Every edge/behavior pair has ties, so MATLAB uses the Student-t rule.
    for threshold_index = 1:no_thresholds
        r_critical = correlation_r_critical( ...
            thresholds(threshold_index), n - 2);
        selected = abs(edge_r) > r_critical;
        pos_by_threshold{threshold_index} = edge_r > 0 & selected;
        neg_by_threshold{threshold_index} = edge_r < 0 & selected;
    end
    return
end

% Completely untied behavior. Most connectome edges are also untied and
% use AS89; any edge with ties uses the Student-t approximation.
edge_has_ties = cache.edge_tieadj > 0;
for threshold_index = 1:no_thresholds
    selected = false(size(edge_r));
    selected(~edge_has_ties) = as89_lookup_selection( ...
        n, D(~edge_has_ties), thresholds(threshold_index));
    if any(edge_has_ties)
        r_critical = correlation_r_critical( ...
            thresholds(threshold_index), n - 2);
        selected(edge_has_ties) = ...
            abs(edge_r(edge_has_ties)) > r_critical;
    end
    pos_by_threshold{threshold_index} = edge_r > 0 & selected;
    neg_by_threshold{threshold_index} = edge_r < 0 & selected;
end
end


function selected = as89_lookup_selection(n, D, threshold)
% Cache the discrete AS89 decision table by n and threshold. Untied rank
% distances are integers, so no per-edge p-value array is needed. This
% branch exists because MATLAB does not use the usual t approximation for a
% completely untied ordinary Spearman pair.
persistent lookup
if isempty(lookup)
    lookup = containers.Map('KeyType', 'char', 'ValueType', 'any');
end
key = sprintf('n%d_t%.17g', n, threshold);
if ~isKey(lookup, key)
    n3const = (n^3 - n) / 3;
    lower_distance = (0:floor(n3const / 2))';
    lookup(key) = spearman_as89_pvalue(n, lower_distance) < threshold;
end
decision = lookup(key);
n3const = (n^3 - n) / 3;
two_sided_distance = min(D, n3const - D);
index = round(two_sided_distance) + 1;
selected = decision(index);
end


function edge_r = cached_linear_correlation(cache, behavior)
if isempty(cache.covars)
    behavior_residual = behavior;
else
    if strcmpi(cache.corr_type, 'Spearman')
        fixed_behavior = tiedrank(behavior);
    else
        fixed_behavior = behavior;
    end
    behavior_residual = fixed_behavior - cache.covar_design * ...
        (cache.covar_design \ fixed_behavior);
    tolerance = max(cache.n, cache.covar_rank) * ...
        eps(class(fixed_behavior)) * ...
        sqrt(sum(abs(fixed_behavior).^2, 1));
    if sqrt(sum(abs(behavior_residual).^2, 1)) < tolerance
        behavior_residual(:) = 0;
    end
end

behavior_centered = behavior_residual - ...
    sum(behavior_residual, 1) / cache.n;
edge_r = (cache.fixed_centered' * behavior_centered) ./ ...
    cache.fixed_norm ./ vecnorm(behavior_centered, 2, 1);
outside = abs(edge_r) > 1;
edge_r(outside) = sign(edge_r(outside));
end


function r_critical = correlation_r_critical(threshold, df)
t_critical = tinv(1 - threshold / 2, df);
r_critical = sqrt(t_critical.^2 ./ (df + t_critical.^2));
end


function pval = spearman_as89_pvalue(n, D)
n3const = (n^3 - n) / 3;
pval = 2 * as89(max(D, n3const - D), n, n3const);
pval(pval > 1) = 1;
end


function p = as89(D, n, n3const)
c = [.2274 .2531 .1745 .0758 .1033 .3932 ...
    .0879 .0151 .0072 .0831 .0131 .00046];
x = (2 * (D - 1) ./ n3const - 1) * sqrt(n - 1);
y = x .* x;
u = x .* (c(1) + (c(2) + c(3)/n)/n + ...
    y .* (-c(4) + (c(5) + c(6)/n)/n - ...
    y .* (c(7) + c(8)/n - y .* (c(9) - c(10)/n + ...
    y .* (c(11) - c(12) * y)/n))/n))/n;
p = u ./ exp(0.5 * y) + 0.5 * erfc(x ./ sqrt(2));
p(p > 1) = 1;
p(p < 0) = 0;
end
