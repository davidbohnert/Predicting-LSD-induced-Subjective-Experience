function [b_comb, b_pos, b_neg] = cpm_fit_three_models_min_norm( ...
        behavior, sum_pos, sum_neg)
% CPM_FIT_THREE_MODELS_MIN_NORM
% Fit the three CPM models with minimum-norm least squares.
%
% LOSO folds can contain an empty positive or negative network. Those
% exact cases are handled explicitly. Clearly full-rank designs use a
% pivoted QR solve; rank-deficient or numerically uncertain designs retain
% PINV and therefore the historical minimum-norm behavior.
% This helper belongs only to the revision LOSO analysis; it does not alter
% the submitted cross-validation or 10-fold consensus model specification.

behavior = behavior(:);
sum_pos = sum_pos(:);
sum_neg = sum_neg(:);
intercept = ones(numel(behavior), 1, 'like', behavior);
pos_empty = ~any(sum_pos);
neg_empty = ~any(sum_neg);

if pos_empty && neg_empty
    intercept_only = sum(behavior) / numel(behavior);
    b_comb = [0; 0; intercept_only];
    b_pos = [0; intercept_only];
    b_neg = [0; intercept_only];
elseif pos_empty
    b_neg = solve_with_min_norm_fallback([sum_neg, intercept], behavior);
    b_comb = [0; b_neg];
    b_pos = [0; sum(behavior) / numel(behavior)];
elseif neg_empty
    b_pos = solve_with_min_norm_fallback([sum_pos, intercept], behavior);
    b_comb = [b_pos(1); 0; b_pos(2)];
    b_neg = [0; sum(behavior) / numel(behavior)];
else
    b_comb = solve_with_min_norm_fallback( ...
        [sum_pos, sum_neg, intercept], behavior);
    b_pos = solve_with_min_norm_fallback([sum_pos, intercept], behavior);
    b_neg = solve_with_min_norm_fallback([sum_neg, intercept], behavior);
end
end


function coefficients = solve_with_min_norm_fallback(design, response)
[Q, R, permutation] = qr(design, 0);
no_columns = size(design, 2);
rank_tolerance = max(size(design)) * eps(R(1));
clearly_full_rank = all(abs(diag(R)) > rank_tolerance);

if clearly_full_rank
    coefficients = zeros(no_columns, 1, 'like', design);
    coefficients(permutation) = R \ (Q' * response);
else
    coefficients = pinv(design) * response;
end
end
