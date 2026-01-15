% DESCRIPTION:
%   Master script for Connectome-based Predictive Modeling (CPM) analysis.
%   Systematically tests parameters (thresholds, k-folds) and matrices
%   (LSD, LSD_difference, LSD_GSR) against behavioral scales.
%
% USAGE:
%   Ensure ./data/ contains:
%       ├─ behav/           (Behavioral .mat vectors)
%       ├─ connectomes/     (Connectivity .mat matrices)
%       └─ covars.mat
%
% ASSUMPTIONS & DATA REQUIREMENTS:
% 1. Subject Sorting:
%    Data is assumed to be pre-sorted. Row i in 'behav/*.mat' MUST correspond
%    to Index i in the 3rd dimension of 'connectomes/*.mat'.
%    The script does NOT match subjects by ID strings.
%
% 2. File Structure (.mat):
%    - Connectivity files must contain a 3D matrix (Node x Node x Subjects).
%    - Behavioral files must contain a single column vector (Subjects x 1).
%    - The script loads the FIRST variable found in each .mat file. Ensure
%      saved .mat files do not contain extraneous variables (like 'ans').
%
% 3. Missing Data:
%    Input matrices and vectors must not contain NaNs or Infs.
%    This pipeline does not perform imputation or row-deletion for missing data.
%
% 4. Fixed Subsets:
%    The variable 'subsetIdx = 20:67' (Line 25) hardcodes the indices for
%    MEQ30 analysis. Ensure your data is sorted such that indices 20:67
%    capture the correct sub-cohort.
%
% 5. Dependencies:
%    - Statistics and Machine Learning Toolbox (for partialcorr, regress)
%    - Parallel Computing Toolbox (for parpool, parfor)
% -------------------------------------------------------------------------

clear; clc;

%% 1. CONFIGURATION & PATH SETUP

scriptDir = fileparts(mfilename('fullpath'));
baseDir = fullfile(scriptDir, '..', '..');

inputPath = fullfile(baseDir, 'data');
outputPath = fullfile(baseDir, 'results');

if ~exist(outputPath, 'dir'), mkdir(outputPath); end

outputTxt = fullfile(outputPath, 'permutation_results.txt');
outputXlsx = fullfile(outputPath, 'permutation_results.xlsx');

% Parameters
BEHAV_SCALES = {'LSD_BDE', 'LSD_GDE', 'LSD_OBN', 'LSD_AED', 'LSD_MEQ30'};
ALT_THRESHOLDS   = {0.05, 0.005};
K_FOLDS      = 10; % Default
BASE_THR     = 0.01; % Default

subsetIdx = 20:67; % Subset of subjects for MEQ30, where there are fewer behavioral ratings

% Parallel Setup
if isempty(gcp('nocreate'))
    parpool();
end

% Load Covariates
covarsFile = fullfile(inputPath, 'covars.mat');
if isfile(covarsFile)
    tmp = load(covarsFile);
    fNames = fieldnames(tmp);
    covarsData = tmp.(fNames{1}); % Assumes structure is inside, or the matrix itself
    if isstruct(covarsData) && isfield(covarsData, 'covars')
        covarsData = covarsData.covars;
    end
else
    error('Covariates file missing: %s', covarsFile);
end

results = {};

%% 2. MAIN ANALYSIS LOOP

for bIdx = 1:length(BEHAV_SCALES)
    behavName = BEHAV_SCALES{bIdx};
    fprintf('\n========================================\n');
    fprintf('Processing Behavior: %s\n', behavName);
    fprintf('========================================\n');

    % --- Load Behavior (Dynamic Name) ---
    behavFile = fullfile(inputPath, 'behav', [behavName '.mat']);
    if ~isfile(behavFile), warning('File not found: %s', behavFile); continue; end

    tmp = load(behavFile);
    fNames = fieldnames(tmp);
    allBehav = tmp.(fNames{1}); % Takes the first variable found

    % --- Load Connectome (Standard LSD) ---
    matrixName = 'LSD';

    covars = covarsData;

    connFile = fullfile(inputPath, 'connectomes', [matrixName '.mat']);
    if ~isfile(connFile), warning('File not found: %s', connFile); continue; end

    tmp = load(connFile);
    fNames = fieldnames(tmp);
    allMats = tmp.(fNames{1});

    % Adjust matrix and covariates for MEQ30, which includes less subjects
    if contains(behavName, 'MEQ30')
        fprintf('   ! Applying subject filter (20:67) for %s...\n', behavName);

        % Slice the Connectome (assuming dim 3 is subjects)
        allMats = allMats(:, :, subsetIdx);

        % Slice the Covariates (assuming dim 1 is subjects)
        covars = covarsData(subsetIdx, :);
    end

    noSub = size(allMats, 3);


    %% BLOCK A: Standard Analysis (k=10, THR=0.01)
    fprintf('-> [Block A] Standard Analysis...\n');
    results = permutation_test_cv(matrixName, behavName, allMats, allBehav, K_FOLDS, BASE_THR, [], outputTxt, results);


    %% BLOCK B: Covariate Analysis
    fprintf('-> [Block B] Covariate Analysis...\n');
    results = permutation_test_cv(matrixName, behavName, allMats, allBehav, K_FOLDS, BASE_THR, covars, outputTxt, results);


    %% BLOCK C: Threshold Sensitivity
    for tIdx = 1:length(ALT_THRESHOLDS)
        currThr = ALT_THRESHOLDS{tIdx};

        fprintf('-> [Block C] Testing Threshold: %.3f\n', currThr);
        % Note: Passing [] for covars as per original logic for sensitivity checks
        results = permutation_test_cv(matrixName, behavName, allMats, allBehav, K_FOLDS, currThr, [], outputTxt, results);
    end


    %% BLOCK D: k-Fold Sensitivity
    kList = [5, noSub];
    for kVal = kList
        fprintf('-> [Block D] Testing k-fold: %d\n', kVal);
        results = permutation_test_cv(matrixName, behavName, allMats, allBehav, kVal, BASE_THR, [], outputTxt, results);
    end


    %% BLOCK E: Difference Matrix Analysis
    diffMatrixName = 'LSD_difference';
    diffConnFile = fullfile(inputPath, 'connectomes', [diffMatrixName '.mat']);

    % Only proceed if the Difference Matrix exists
    if isfile(diffConnFile)

        % --- 1. Determine which Behavior Vector to use ---
        % Construct the expected filename for difference behavior
        targetDiffBehavName = [behavName '_difference'];
        diffBehavFile = fullfile(inputPath, 'behav', [targetDiffBehavName '.mat']);

        if isfile(diffBehavFile)
            tmpB = load(diffBehavFile);
            fNamesB = fieldnames(tmpB);
            diffBehav = tmpB.(fNamesB{1});
            useBehavName = targetDiffBehavName; % Use the specific name for results
        else
            % if not found -> Fallback to standard behavior
            diffBehav = allBehav;
            useBehavName = behavName; % Use the standard name for results
        end

        % --- 2. Load Difference Connectome ---
        tmpC = load(diffConnFile);
        fNamesC = fieldnames(tmpC);
        diffMats = tmpC.(fNamesC{1});

        % Apply Subject Filter if needed
        if contains(behavName, 'MEQ30')
            diffMats = diffMats(:, :, subsetIdx);
        end

        fprintf('-> [Block E] Difference Matrix Analysis...\n');

        % Run Analysis
        results = permutation_test_cv(diffMatrixName, useBehavName, diffMats, diffBehav, K_FOLDS, BASE_THR, [], outputTxt, results);
        results = permutation_test_cv(diffMatrixName, useBehavName, diffMats, diffBehav, K_FOLDS, BASE_THR, covars, outputTxt, results);
    end



    %% BLOCK F: Global Signal Regression Check (LSD_GSR)
    gsrMatrixName = 'LSD_GSR';
    gsrConnFile = fullfile(inputPath, 'connectomes', [gsrMatrixName '.mat']);

    if isfile(gsrConnFile)
        fprintf('-> [Block F] Checking GSR Control (LSD_GSR)...\n');

        tmpGSR = load(gsrConnFile);
        fNamesGSR = fieldnames(tmpGSR);
        gsrMats = tmpGSR.(fNamesGSR{1});

        if contains(behavName, 'MEQ30')
            gsrMats = gsrMats(:, :, subsetIdx);
        end

        results = permutation_test_cv(gsrMatrixName, behavName, gsrMats, allBehav, K_FOLDS, BASE_THR, [], outputTxt, results);
    else
        fprintf('-> [Block F] Skipped (LSD_GSR.mat not found).\n');
    end

end

%% 3. SAVE FINAL RESULTS
fprintf('\nPipeline complete. Saving results...\n');
header = {'Analysis', 'k', 'thresh', 'r_comb', 'p_comb', 'r_pos', 'p_pos', 'r_neg', 'p_neg', ...
    'mse_comb', 'mse_p_comb', 'mse_pos', 'mse_p_pos', 'mse_neg', 'mse_p_neg'};
resultsWithHeader = [header; results];
writecell(resultsWithHeader, outputXlsx);