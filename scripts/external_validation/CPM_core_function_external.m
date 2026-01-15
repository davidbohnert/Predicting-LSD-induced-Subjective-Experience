function [R_comb, R_pos, R_neg, mse_comb, mse_pos, mse_neg, final_pos_mask, final_neg_mask] = ...
    CPM_core_function_external(mats_train, behav_train, ...
    mats_test, behav_test, fixed_partition, covars, corr_type)
% -------------------------------------------------------------------------
% CPM_EXTERNAL_CORE_FUNCTION (formerly predict_external_robust)
%   Trains a CPM model on a Training Dataset using a "Consensus Mask" 
%   approach (edges present in k-folds), then tests it on an independent 
%   External Dataset.
%
% INPUTS:
%   mats_train      - Training Connectivity matrices (Nodes x Nodes x Subj)
%   behav_train     - Training Behavior scores
%   mats_test       - External Test Connectivity matrices
%   behav_test      - External Test Behavior scores
%   covars          - Covariates for Training data (can be empty)
%   fixed_partition - Fixed Partition for Cross-Validation
%   corr_type       - 'Pearson' or 'Spearman'
%
% OUTPUTS:
%   R_comb, etc.    - Correlations on the External Test Set
% -------------------------------------------------------------------------

    %% 1. CONFIGURATION
    P_THR   = 0.01;
    CONSENSUS_THR = 0.8; % Edge must appear in 80% of folds to be kept
    
    
    %% 2. INITIALIZATION
    [no_node, ~, num_train] = size(mats_train);
    num_test = size(mats_test, 3);
    
    % Prepare Consensus Masks
    sum_pos_mask = zeros(no_node, no_node);
    sum_neg_mask = zeros(no_node, no_node);
    
    % Create CV Partition for Consensus Building
    % Note: RNG is controlled by the external wrapper script
    k = fixed_partition.NumTestSets;


    %% 3. CONSENSUS MASK GENERATION (Training Set)
    
    for fold = 1:k
        
        % --- Define Fold ---
        train_idx = training(fixed_partition, fold);
        
        % Extract Data for this fold
        curr_mats  = mats_train(:, :, train_idx);
        curr_behav = behav_train(train_idx);
        
        % Flatten for Correlation
        curr_vcts  = reshape(curr_mats, [], sum(train_idx));
        curr_behav = curr_behav(:);
        
        % --- Feature Selection ---
        if isempty(covars)
            [r_mat, p_mat] = corr(curr_vcts', curr_behav, 'type', corr_type);
        else
            curr_covars = covars(train_idx, :);
            [r_mat, p_mat] = partialcorr(curr_vcts', curr_behav, curr_covars, 'type', corr_type);
        end
        
        % Reshape to Matrix
        r_mat = reshape(r_mat, no_node, no_node);
        p_mat = reshape(p_mat, no_node, no_node);
        
        % Define Edges
        pos_edges = (r_mat > 0 & p_mat < P_THR);
        neg_edges = (r_mat < 0 & p_mat < P_THR);
        
        % Accumulate (Vote)
        sum_pos_mask = sum_pos_mask + pos_edges;
        sum_neg_mask = sum_neg_mask + neg_edges;
    end
    
    % Final Consensus: Keep edges that appeared in >80% of folds
    final_pos_mask = sum_pos_mask >= (CONSENSUS_THR * k);
    final_neg_mask = sum_neg_mask >= (CONSENSUS_THR * k);

    %% 4. MODEL TRAINING (Full Training Set)
    
    % Calculate Sums for ALL Training Subjects using Final Mask
    % Optimization: Vectorized Dot Product
    % (Mask Vector * All_Subjects_Matrix) / 2
    
    train_vcts_all = reshape(mats_train, [], num_train);
    
    train_sumpos = (final_pos_mask(:)' * train_vcts_all)' / 2;
    train_sumneg = (final_neg_mask(:)' * train_vcts_all)' / 2;
    
    % Regress Behavior on Sums
    % Suppress rank warning (common in CPM if sums are collinear)
    warning('off', 'stats:regress:RankDefDesignMat');
    
    behav_train = behav_train(:);
    
    % Model 1: Combined
    b_comb = regress(behav_train, [train_sumpos, train_sumneg, ones(num_train, 1)]);
    
    % Model 2: Positive
    b_pos  = regress(behav_train, [train_sumpos, ones(num_train, 1)]);
    
    % Model 3: Negative
    b_neg  = regress(behav_train, [train_sumneg, ones(num_train, 1)]);
    
    warning('on', 'stats:regress:RankDefDesignMat');


    %% 5. EXTERNAL PREDICTION (Test Set)
    
    % Calculate Sums for External Test Subjects
    test_vcts_all = reshape(mats_test, [], num_test);
    
    test_sumpos = (final_pos_mask(:)' * test_vcts_all)' / 2;
    test_sumneg = (final_neg_mask(:)' * test_vcts_all)' / 2;
    
    % Apply Training Weights (b) to Test Data
    behav_pred_comb = b_comb(1)*test_sumpos + b_comb(2)*test_sumneg + b_comb(3);
    behav_pred_pos  = b_pos(1)*test_sumpos  + b_pos(2);
    behav_pred_neg  = b_neg(1)*test_sumneg  + b_neg(2);


    %% 6. EVALUATION
    
    behav_test = behav_test(:);
    
    % Correlations
    [R_comb, ~] = corr(behav_pred_comb, behav_test);
    [R_pos, ~]  = corr(behav_pred_pos,  behav_test);
    [R_neg, ~]  = corr(behav_pred_neg,  behav_test);
    
    % MSE
    mse_comb = mean((behav_test - behav_pred_comb).^2);
    mse_pos  = mean((behav_test - behav_pred_pos).^2);
    mse_neg  = mean((behav_test - behav_pred_neg).^2);

end