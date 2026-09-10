function cache = cpm_prepare_correlation_cache(edge_data, covars, corr_type)
% CPM_PREPARE_CORRELATION_CACHE
% Prepare the fixed connectome side of one CPM training fold.
%
% Edge selection remains Pearson, partial Pearson, Spearman, or partial
% Spearman as requested. Fixed edge ranks, residuals, centering, and norms
% are calculated once rather than once per behavior permutation.
%
% In the submitted core, CORR/PARTIALCORR recomputed both sides for every
% shuffle. Within a fixed fold the connectome values and covariates never
% change; only behavior is permuted. This cache stores only those invariant
% quantities and leaves behavior-specific correlation to the selection step.

cache.corr_type = char(validatestring(corr_type, ...
    {'Pearson', 'Spearman'}));
cache.covars = covars;
cache.n = size(edge_data, 1);
cache.covar_design = [];
cache.covar_rank = 0;
cache.edge_data = [];
cache.edge_ranks = [];
cache.edge_tieadj = [];
cache.edge_rank_sumsq = [];
cache.fixed_centered = [];
cache.fixed_norm = [];

is_spearman = strcmpi(cache.corr_type, 'Spearman');
if is_spearman
    [fixed_values, cache.edge_tieadj] = tiedrank(edge_data);
else
    fixed_values = edge_data;
end

if isempty(covars)
    if is_spearman
        % The D statistic gives MATLAB-equivalent Spearman coefficients and
        % preserves its AS89 rule for the occasional completely untied fold.
        cache.edge_ranks = fixed_values;
        cache.edge_tieadj = cache.edge_tieadj(:);
        cache.edge_rank_sumsq = sum(fixed_values.^2, 1)';
    else
        cache.fixed_centered = fixed_values - ...
            sum(fixed_values, 1) / cache.n;
        cache.fixed_norm = vecnorm(cache.fixed_centered, 2, 1)';
    end
else
    if is_spearman
        fixed_covars = tiedrank(covars);
    else
        fixed_covars = covars;
    end

    cache.covar_rank = rank(fixed_covars);
    cache.covar_design = [ones(cache.n, 1), fixed_covars];
    % Partial correlation is ordinary correlation between residuals after
    % removing the same covariate design from edges and behavior.
    fixed_residuals = fixed_values - cache.covar_design * ...
        (cache.covar_design \ fixed_values);

    % Match partialcorr's handling of a fixed edge that is effectively
    % constant after the covariates have been removed.
    tolerance = max(cache.n, cache.covar_rank) * ...
        eps(class(fixed_values)) .* ...
        sqrt(sum(abs(fixed_values).^2, 1));
    residual_norm = sqrt(sum(abs(fixed_residuals).^2, 1));
    fixed_residuals(:, residual_norm < tolerance) = 0;

    cache.fixed_centered = fixed_residuals - ...
        sum(fixed_residuals, 1) / cache.n;
    cache.fixed_norm = vecnorm(cache.fixed_centered, 2, 1)';
end

% Only the very-small-n Spearman fallback needs the original edge values.
% No paper analysis enters this branch.
if is_spearman && cache.n < 10
    cache.edge_data = edge_data;
end
end
