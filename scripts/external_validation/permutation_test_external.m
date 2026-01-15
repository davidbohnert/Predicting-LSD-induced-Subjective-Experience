function results = permutation_test_external(matrix_train_name, behav_train_name, mats_train, behav_train, matrix_test_name, mats_test, behav_test, covars, output_txt, results)
% -------------------------------------------------------------------------
% RUN_EXTERNAL_PERMUTATION
%   Runs permutation testing for External Validation.
%   Trains a CPM model on the training set (e.g., LSD) and tests it on
%   a completely independent dataset (e.g., Other Psych).
%   Significance is determined by shuffling the *training* behavior labels.
%
% INPUTS:
%   matrix_train_name - Name of training dataset (logging)
%   mats_train        - Training connectivity matrices
%   behav_train       - Training behavior scores
%   matrix_test_name  - Name of testing dataset (logging)
%   mats_test         - Testing connectivity matrices
%   behav_test        - Testing behavior scores
%   covars            - Covariates for TRAINING set (can be empty)
%
% OUTPUTS:
%   results           - Updated cell array with stats
% -------------------------------------------------------------------------

%% 1. CONFIGURATION
NO_ITERATIONS = 1000;

k = 10;

rng(123);
no_sub = size(mats_train, 3);

% Creating reproducible partition for consensus feature selection
fixed_partition = cvpartition(no_sub, 'KFold', k);

% --- Covariates & Logging Setup ---
covar_text = '';
if ~isempty(covars)
    covar_text = ' - with covariates';
end

fprintf('Running External Permutation: %s vs %s%s -> Test: %s\n', ...
    matrix_train_name, behav_train_name, covar_text, matrix_test_name);

% --- Determine Correlation Type ---
if contains(behav_train_name, 'LSD_OBN') || contains(behav_train_name, 'LSD_MEQ30')
    corr_type = 'Pearson';
else
    corr_type = 'Spearman';
end


%% 2. TRUE MODEL PREDICTION (Iteration 1)

[t_r_comb, t_r_pos, t_r_neg, t_mse_comb, t_mse_pos, t_mse_neg, final_pos_mask, final_neg_mask] = ...
    CPM_core_function_external(mats_train, behav_train, mats_test, behav_test, fixed_partition, covars, corr_type);

true_r   = [t_r_comb, t_r_pos, t_r_neg];
true_mse = [t_mse_comb, t_mse_pos, t_mse_neg];

% Optional: Export consensus masks for later visualization
edgesDir = fullfile(fileparts(output_txt), 'edges');
if ~exist(edgesDir, 'dir'), mkdir(edgesDir); end

posFile = fullfile(edgesDir, [behav_train_name '_pos.csv']);
negFile = fullfile(edgesDir, [behav_train_name '_neg.csv']);

% Export only if not already present
if ~isfile(posFile), writematrix(final_pos_mask, posFile); end
if ~isfile(negFile), writematrix(final_neg_mask, negFile); end

%% 3. PERMUTATION LOOP

% Pre-allocate full size (including row 1 for True value)
perm_r   = zeros(NO_ITERATIONS, 3);
perm_mse = zeros(NO_ITERATIONS, 3);

% Assign True values to first row
perm_r(1, :)   = true_r;
perm_mse(1, :) = true_mse;

parfor it = 2:NO_ITERATIONS
    % Feedback every 100 iterations to avoid cluttering console
    if mod(it, 100) == 0
        fprintf('  ...Iteration %d / %d\n', it, NO_ITERATIONS);
    end

    % Shuffle TRAINING behavior only
    rng(it + 200);
    shuffled_train_behav = behav_train(randperm(length(behav_train)));

    % Run CPM on Null Training Data -> Predict Real Test Data
    [p_r_comb, p_r_pos, p_r_neg, p_mse_comb, p_mse_pos, p_mse_neg] = ...
        CPM_core_function_external(mats_train, shuffled_train_behav, mats_test, behav_test, fixed_partition, covars, corr_type);

    % Store Null Results
    perm_r(it, :)   = [p_r_comb, p_r_pos, p_r_neg];
    perm_mse(it, :) = [p_mse_comb, p_mse_pos, p_mse_neg];
end


%% 4. CALCULATE P-VALUES

% Correlation (One-tailed: Count if Null >= True)
pval_comb = sum(perm_r(:,1) >= t_r_comb) / NO_ITERATIONS;
pval_pos  = sum(perm_r(:,2) >= t_r_pos)  / NO_ITERATIONS;
pval_neg  = sum(perm_r(:,3) >= t_r_neg)  / NO_ITERATIONS;

% MSE (One-tailed: Count if Null <= True)
pval_mse_comb = sum(perm_mse(:,1) <= t_mse_comb) / NO_ITERATIONS;
pval_mse_pos  = sum(perm_mse(:,2) <= t_mse_pos)  / NO_ITERATIONS;
pval_mse_neg  = sum(perm_mse(:,3) <= t_mse_neg)  / NO_ITERATIONS;


%% 5. SAVE & RETURN RESULTS

% Log to Text File
try
    fileID = fopen(output_txt, 'a');
    fprintf(fileID, '\n------------------------------------------------\n');
    fprintf(fileID, '%s vs %s%s\nTest Set: %s\n', matrix_train_name, behav_train_name, covar_text, matrix_test_name);

    fprintf(fileID, 'R_comb = %1.4f (p = %1.4f)\n', t_r_comb, pval_comb);
    fprintf(fileID, 'R_pos  = %1.4f (p = %1.4f)\n', t_r_pos,  pval_pos);
    fprintf(fileID, 'R_neg  = %1.4f (p = %1.4f)\n', t_r_neg,  pval_neg);

    fprintf(fileID, 'MSE_comb = %1.4f (p = %1.4f)\n', t_mse_comb, pval_mse_comb);
    fprintf(fileID, 'MSE_pos  = %1.4f (p = %1.4f)\n', t_mse_pos,  pval_mse_pos);
    fprintf(fileID, 'MSE_neg  = %1.4f (p = %1.4f)\n', t_mse_neg,  pval_mse_neg);
    fclose(fileID);
catch ME
    warning('CPM:LogWriteFailed', 'Failed to write to text file: %s', ME.message);
end

% Append to Results Cell Array
results(end+1, :) = {sprintf('%s vs %s%s', matrix_train_name, behav_train_name, covar_text), matrix_test_name, ...
    t_r_comb, pval_comb, t_r_pos, pval_pos, t_r_neg, pval_neg, ...
    t_mse_comb, pval_mse_comb, t_mse_pos, pval_mse_pos, t_mse_neg, pval_mse_neg};

end
