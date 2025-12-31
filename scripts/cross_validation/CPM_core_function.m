function [R_comb, R_pos, R_neg, mse_comb, mse_pos, mse_neg] = CPM_core_function(all_mats, all_behav, cv_partition, THR, covars, corr_type)
% -------------------------------------------------------------------------
% CPM_CORE_FUNCTION
%   Performs Connectome-based Predictive Modeling (CPM) using a fixed
%   cross-validation partition.
%
% INPUTS:
%   all_mats     - Connectivity matrices (Nodes x Nodes x Subjects)
%   all_behav    - Behavioral scores (Subjects x 1)
%   cv_partition - A cvpartition object defining the folds
%   THR          - P-value threshold for edge selection (e.g., 0.01)
%   covars       - Covariates matrix (Subjects x Regressors), or []
%   corr_type    - 'Pearson' or 'Spearman'
%
% OUTPUTS:
%   R_comb, etc. - Correlation between predicted and observed scores
%   mse_comb, etc.- Mean Squared Error of predictions
% -------------------------------------------------------------------------

    %% 1. INITIALIZATION
    no_sub = size(all_mats, 3);
    no_node = size(all_mats, 1);
    
    % Derive k (number of folds) from the partition object
    k = cv_partition.NumTestSets;

    % Pre-allocate prediction vectors
    behav_pred_pos = zeros(no_sub, 1);
    behav_pred_neg = zeros(no_sub, 1);
    behav_pred_comb = zeros(no_sub, 1); % Corrected name to match output

    %% 2. CROSS-VALIDATION LOOP
    for fold = 1:k
        
        % --- Define Train/Test Sets ---
        % CRITICAL FIX: Wrap in find() to convert logical mask (0/1) to indices (1, 2, 5...)
        train_idx = find(training(cv_partition, fold));
        test_idx  = find(test(cv_partition, fold));

        % Extract Training Data
        train_mats = all_mats(:,:,train_idx);
        train_behav = all_behav(train_idx);
        
        % Flatten matrices: (Nodes x Nodes x Subj) -> (Edges x Subj)
        train_vcts = reshape(train_mats, [], length(train_idx));
        
        % --- Feature Selection (Correlation) ---
        train_behav = train_behav(:); % Ensure column vector
        
        edge_no = size(train_vcts, 1);
        
        if isempty(covars)
            [r_mat, p_mat] = corr(train_vcts', train_behav, 'type', corr_type);
        else
            train_covars = covars(train_idx, :);
            [r_mat, p_mat] = partialcorr(train_vcts', train_behav, train_covars, 'type', corr_type);
        end

        % Reshape vectors back to matrices for mask creation
        r_mat = reshape(r_mat, no_node, no_node);
        p_mat = reshape(p_mat, no_node, no_node);

        % Create Binary Masks
        pos_mask = zeros(no_node, no_node);
        neg_mask = zeros(no_node, no_node);

        pos_edge = find(r_mat > 0 & p_mat < THR);
        neg_edge = find(r_mat < 0 & p_mat < THR);

        pos_mask(pos_edge) = 1;
        neg_mask(neg_edge) = 1;

        % --- Model Building (Training Set) ---
        num_train = length(train_idx);
        train_sumpos = zeros(num_train, 1);
        train_sumneg = zeros(num_train, 1);

        for ss = 1:num_train
            % Sum edges for each training subject (divided by 2 for symmetry)
            train_sumpos(ss) = sum(sum(train_mats(:,:,ss) .* pos_mask)) / 2;
            train_sumneg(ss) = sum(sum(train_mats(:,:,ss) .* neg_mask)) / 2;
        end

        % Fit Linear Regressions
        % Model 1: Combined (Pos + Neg + Intercept)
        b_comb = regress(train_behav, [train_sumpos, train_sumneg, ones(num_train, 1)]);
        
        % Model 2: Positive Only (Pos + Intercept)
        b_pos  = regress(train_behav, [train_sumpos, ones(num_train, 1)]);
        
        % Model 3: Negative Only (Neg + Intercept)
        b_neg  = regress(train_behav, [train_sumneg, ones(num_train, 1)]);

        
        % --- Model Testing (Test Set) ---
        test_mats = all_mats(:,:,test_idx);
        num_test = length(test_idx);

        for i = 1:num_test
            % Calculate strength for test subject
            test_sumpos = sum(sum(test_mats(:,:,i) .* pos_mask)) / 2;
            test_sumneg = sum(sum(test_mats(:,:,i) .* neg_mask)) / 2;

            % Apply the training model parameters to the test subject
            % Note: test_idx(i) ensures we map prediction back to the correct subject ID
            behav_pred_pos(test_idx(i)) = b_pos(1)*test_sumpos + b_pos(2);
            behav_pred_neg(test_idx(i)) = b_neg(1)*test_sumneg + b_neg(2);
            
            % Combined prediction
            behav_pred_comb(test_idx(i)) = b_comb(1)*test_sumpos + b_comb(2)*test_sumneg + b_comb(3);
        end
    end

    %% 3. COMPARE PREDICTED VS OBSERVED
    all_behav = all_behav(:);
    
    [R_pos, ~]  = corr(behav_pred_pos, all_behav);
    [R_neg, ~]  = corr(behav_pred_neg, all_behav);
    [R_comb, ~] = corr(behav_pred_comb, all_behav);

    mse_comb = mean((all_behav - behav_pred_comb).^2);
    mse_pos  = mean((all_behav - behav_pred_pos).^2);
    mse_neg  = mean((all_behav - behav_pred_neg).^2);
    
end