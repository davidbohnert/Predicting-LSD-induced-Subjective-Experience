function results = permutation_test_cv(matrix_name, behav_name, all_mats, all_behav, k, THR, covars, output_txt, results)
% -------------------------------------------------------------------------
% RUN_PERMUTATION_TEST
%   Runs 1000 permutation tests to validate CPM models. 
%   Compares True Prediction (r, MSE) against a null distribution generated 
%   by shuffling behavioral labels.
%
% INPUTS:
%   matrix_name  - String name of the connectivity matrix (for logging)
%   behav_name   - String name of the behavior (for logging)
%   all_mats     - Connectivity matrices (Subjects x Nodes x Nodes or Vectorized)
%   all_behav    - Behavioral scores (Subjects x 1)
%   k            - Number of folds for Cross-Validation
%   THR          - P-value threshold for feature selection
%   covars       - Covariates matrix (can be empty)
%   output_txt   - Path to log text file
%   results      - Cell array to append final results to
%
% OUTPUTS:
%   results      - Updated cell array with new statistics
% -------------------------------------------------------------------------

    %% 1. SETUP & CONFIGURATION
    NO_ITERATIONS = 1000; % Standard for permutation testing

    % Setting the seed for cross validation for reproducibale folds
    rng(123); 
    no_sub = size(all_mats, 3);
    
    cv_fixed = cvpartition(no_sub, 'KFold', k);
    
    % --- Handle Covariates & Special Cases ---
    covar_text = '';
    if ~isempty(covars)
        covar_text = ' - with covariates';
    end

    % --- Determine Correlation Type ---
    % Use Pearson for OBN/MEQ30, Spearman for others (based on data distribution)
    if contains(behav_name, 'LSD_OBN') || contains(behav_name, 'LSD_MEQ30')
        corr_type = 'Pearson';
    else
        corr_type = 'Spearman';
    end

    fprintf('Running Permutation: %s vs %s%s (k=%d, p=%.3f)\n', ...
        matrix_name, behav_name, covar_text, k, THR);


    %% 2. TRUE MODEL PREDICTION (Iteration 1)
    
    [t_r_comb, t_r_pos, t_r_neg, t_mse_comb, t_mse_pos, t_mse_neg] = ...
        CPM_core_function(all_mats, all_behav, cv_fixed, THR, covars, corr_type);

    % Store True Results
    true_r   = [t_r_comb, t_r_pos, t_r_neg];
    true_mse = [t_mse_comb, t_mse_pos, t_mse_neg];


    %% 3. PERMUTATION LOOP (Iterations 2 -> NO_ITERATIONS)
    
    % Pre-allocate full matrices (include row 1 for True value)
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
        
        % Shuffle Behavior (Null Model)
        % rng(it) ensures reproducibility across parallel workers
        rng(it + 200);
        shuffled_behav = all_behav(randperm(length(all_behav)));

        % Run CPM on Null Data
        [p_r_comb, p_r_pos, p_r_neg, p_mse_comb, p_mse_pos, p_mse_neg] = ...
            CPM_core_function(all_mats, shuffled_behav, cv_fixed, THR, covars, corr_type);

        % Store Null Results
        perm_r(it, :)   = [p_r_comb, p_r_pos, p_r_neg];
        perm_mse(it, :) = [p_mse_comb, p_mse_pos, p_mse_neg];
    end


    %% 4. CALCULATE P-VALUES
    % (Count how many null results are "better" than the true result)
    
    % For R-values (Higher is better, define p as count >= true)
    pval_comb = sum(perm_r(:,1) >= t_r_comb) / NO_ITERATIONS;
    pval_pos  = sum(perm_r(:,2) >= t_r_pos)  / NO_ITERATIONS;
    pval_neg  = sum(perm_r(:,3) >= t_r_neg)  / NO_ITERATIONS;
    
    % For MSE (Lower is better, define p as count <= true)
    pval_mse_comb = sum(perm_mse(:,1) <= t_mse_comb) / NO_ITERATIONS;
    pval_mse_pos  = sum(perm_mse(:,2) <= t_mse_pos)  / NO_ITERATIONS;
    pval_mse_neg  = sum(perm_mse(:,3) <= t_mse_neg)  / NO_ITERATIONS;


    %% 5. SAVE & RETURN RESULTS
    
    % Append to text file log
    try
        fileID = fopen(output_txt, 'a');
        fprintf(fileID, '\n------------------------------------------------\n');
        fprintf(fileID, '%s vs %s%s\nParams: k=%d, THR=%.3f\n', matrix_name, behav_name, covar_text, k, THR);
        
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
    results(end+1, :) = {sprintf('%s vs %s%s', matrix_name, behav_name, covar_text), ...
               k, THR, ...
               t_r_comb, pval_comb, t_r_pos, pval_pos, t_r_neg, pval_neg, ...
               t_mse_comb, pval_mse_comb, t_mse_pos, pval_mse_pos, t_mse_neg, pval_mse_neg};

end
