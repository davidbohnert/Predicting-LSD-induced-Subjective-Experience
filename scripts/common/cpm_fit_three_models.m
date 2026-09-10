function [b_comb, b_pos, b_neg] = cpm_fit_three_models( ...
        behavior, sum_pos, sum_neg)
% CPM_FIT_THREE_MODELS
% Fit the conventional combined, positive, and negative CPM regressions.
%
% Only coefficients are needed in the permutation pipelines. The local
% solver follows REGRESS's coefficient-only, rank-revealing QR calculation
% without invoking its repeated argument, missing-data, and warning setup.
% Upstream CPM input validation guarantees finite, aligned vectors.
% The three design matrices and coefficient order are identical to the
% combined, positive-only, and negative-only REGRESS calls in the submitted
% core functions.

behavior = behavior(:);
sum_pos = sum_pos(:);
sum_neg = sum_neg(:);
intercept = ones(numel(behavior), 1, 'like', behavior);

b_comb = coefficient_only_regress( ...
    behavior, [sum_pos, sum_neg, intercept]);
b_pos = coefficient_only_regress(behavior, [sum_pos, intercept]);
b_neg = coefficient_only_regress(behavior, [sum_neg, intercept]);
end


function coefficients = coefficient_only_regress(response, design)
% Match the one-output path used by MATLAB REGRESS for finite CPM inputs,
% including zero coefficients for columns rejected as rank-dependent.
[no_rows, no_columns] = size(design);
[Q, R, permutation] = qr(design, 0);

if isempty(R)
    design_rank = 0;
elseif isvector(R)
    design_rank = double(abs(R(1)) > 0);
else
    design_rank = sum(abs(diag(R)) > ...
        max(no_rows, no_columns) * eps(R(1)));
end

if design_rank < no_columns
    R = R(1:design_rank, 1:design_rank);
    Q = Q(:, 1:design_rank);
    permutation = permutation(1:design_rank);
end

prototype = design(1) + response(1);
coefficients = zeros(no_columns, 1, 'like', prototype);
coefficients(permutation) = R \ (Q' * response);
end
