% DESCRIPTION:
%   Master script for External Validation of CPM models.
%   Trains models on a primary dataset (e.g., LSD) and tests their predictive 
%   power on independent external datasets (e.g., 'other_psych', 'all_amphs').
%   Utilizes a "Consensus Mask" approach where edges must appear in k-folds
%   of the training set to be selected.
%
% USAGE:
%   Ensure ./data/ contains:
%       ├─ behav/           (Behavioral .mat vectors for both Train & Test)
%       ├─ connectomes/     (Connectivity .mat matrices for both Train & Test)
%       └─ covars.mat       (Covariates for the Training set)
%
% ASSUMPTIONS & DATA REQUIREMENTS:
% 1. Independent Subject Sorting: 
%    Data is assumed to be pre-sorted. 
%    - Training Set: Row i in 'LSD_*.mat' must match Index i in 'LSD.mat'.
%    - Testing Set:  Row j in 'Test_*.mat' must match Index j in 'Test.mat'.
%    The script does NOT match subjects by ID strings across files.
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
% 4. Fixed Subsets (MEQ30):
%    The variable 'subsetIdx = 20:67' (Configuration) hardcodes the indices 
%    for MEQ30 analysis. This assumes the TRAINING data (LSD) is sorted 
%    such that indices 20:67 capture the subjects who have MEQ30 scores.
%
% 5. Dependencies:
%    - Statistics and Machine Learning Toolbox (for partialcorr, regress)
%    - Parallel Computing Toolbox (for parpool, parfor)
% -------------------------------------------------------------------------

clear; clc;

%% 1. CONFIGURATION & PATH SETUP

scriptDir = fileparts(mfilename('fullpath'));
baseDir = fullfile(scriptDir, '..'); 

inputPath = fullfile(baseDir, 'data');
outputPath = fullfile(baseDir, 'results');

if ~exist(outputPath, 'dir'), mkdir(outputPath); end

outputTxt = fullfile(outputPath, 'external_permutation_results.txt');
outputXlsx = fullfile(outputPath, 'external_permutation_results.xlsx');

% --- Analysis Parameters ---
% Define exactly which behaviors and datasets to process here.
BEHAV_SCALES = {'good', 'AED', 'MEQ30'};      % e.g., {'good', 'AED', 'MEQ30'}
TEST_DATASETS = {'other_psych', 'all_amphs'}; % e.g., {'other_psych', 'all_amphs'}

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
    covarsData = tmp.(fNames{1}); % Assumes structure or matrix
    if isstruct(covarsData) && isfield(covarsData, 'covars')
         covarsData = covarsData.covars;
    end
else
    error('Covariates file missing: %s', covarsFile);
end

results = {};

%% 2. MAIN VALIDATION LOOP

for tIdx = 1:length(TEST_DATASETS)
    testDataName = TEST_DATASETS{tIdx};
    
    for bIdx = 1:length(BEHAV_SCALES)
        behavSuffix = BEHAV_SCALES{bIdx};
        
        fprintf('\n========================================\n');
        fprintf('External Validation: %s -> %s\n', testDataName, behavSuffix);
        fprintf('========================================\n');
        
        %% PREPARE DATA
        
        % 1. Construct File Names
        matrixTrainName = 'LSD';
        behavTrainName  = ['LSD_' behavSuffix];
        behavTestName   = [testDataName '_' behavSuffix];
        
        % 2. Load Behavior Vectors
        % Train
        fileB_Tr = fullfile(inputPath, 'behav', [behavTrainName '.mat']);
        if ~isfile(fileB_Tr), warning('Missing Train Behav: %s', fileB_Tr); continue; end
        tmp = load(fileB_Tr); f = fieldnames(tmp);
        behavTrain = tmp.(f{1});
        
        % Test
        fileB_Te = fullfile(inputPath, 'behav', [behavTestName '.mat']);
        if ~isfile(fileB_Te), warning('Missing Test Behav: %s', fileB_Te); continue; end
        tmp = load(fileB_Te); f = fieldnames(tmp);
        behavTest = tmp.(f{1});
        
        covars = covarsData;

        matrixTestName = testDataName;
        
        % 4. Load Connectivity Matrices
        % Train
        fileM_Tr = fullfile(inputPath, 'connectomes', [matrixTrainName '.mat']);
        if ~isfile(fileM_Tr), warning('Missing Train Matrix: %s', fileM_Tr); continue; end
        tmp = load(fileM_Tr); f = fieldnames(tmp);
        matsTrain = tmp.(f{1});
        
        % Test
        fileM_Te = fullfile(inputPath, 'connectomes', [matrixTestName '.mat']);
        if ~isfile(fileM_Te), warning('Missing Test Matrix: %s', fileM_Te); continue; end
        tmp = load(fileM_Te); f = fieldnames(tmp);
        matsTest = tmp.(f{1});

        % Adjust matrix and covariates for MEQ30, which includes less subjects
        if contains(behavTrainName, 'MEQ30')
            fprintf('   ! Applying subject filter (20:67) for %s...\n', behavTrainName);

            % Slice the Connectome (assuming dim 3 is subjects)
            matsTrain = matsTrain(:, :, subsetIdx);

            % Slice the Covariates (assuming dim 1 is subjects)
            covars = covarsData(subsetIdx, :);
        end
        
        %% BLOCK A: Standard Analysis (No Covariates)
        fprintf('-> [Block A] Standard Analysis...\n');
        results = permutation_test_external(matrixTrainName, behavTrainName, matsTrain, behavTrain, ...
                                           matrixTestName, matsTest, behavTest, [], outputTxt, results);
                                       
        

        %% BLOCK B: Covariate Analysis
        fprintf('-> [Block B] Covariate Analysis...\n');
        results = permutation_test_external(matrixTrainName, behavTrainName, matsTrain, behavTrain, ...
                                           matrixTestName, matsTest, behavTest, covars, outputTxt, results);
        
        
        %% BLOCK C: Difference Matrix Analysis
        fprintf('-> [Block C] Difference Matrix Analysis...\n');
        
        diffTrainMatrix = 'LSD_difference';
        
        % Construct Difference File Names
        diffTestMatrix = [matrixTestName '_difference'];
        diffTrainBehav = [behavTrainName '_difference'];
        diffTestBehav  = [behavTestName '_difference'];
        
        % Check Files Existance
        pathM_Tr = fullfile(inputPath, 'connectomes', [diffTrainMatrix '.mat']);
        pathM_Te = fullfile(inputPath, 'connectomes', [diffTestMatrix '.mat']);
        pathB_Tr = fullfile(inputPath, 'behav', [diffTrainBehav '.mat']);
        pathB_Te = fullfile(inputPath, 'behav', [diffTestBehav '.mat']);
        
        if isfile(pathM_Tr) && isfile(pathM_Te) && isfile(pathB_Tr) && isfile(pathB_Te)
            
            % Load Difference Data
            tmp = load(pathM_Tr); f=fieldnames(tmp); dMatsTr = tmp.(f{1});
            tmp = load(pathM_Te); f=fieldnames(tmp); dMatsTe = tmp.(f{1});
            tmp = load(pathB_Tr); f=fieldnames(tmp); dBehavTr = tmp.(f{1});
            tmp = load(pathB_Te); f=fieldnames(tmp); dBehavTe = tmp.(f{1});


            % Adjust for MEQ30
            if contains(behavTrainName, 'LSD_MEQ30')
                dMatsTr = dMatsTr(:, :, subsetIdx);
            end
            
            % 1. Without Covariates
            results = permutation_test_external(diffTrainMatrix, diffTrainBehav, dMatsTr, dBehavTr, ...
                                               diffTestMatrix, dMatsTe, dBehavTe, [], outputTxt, results);
            
            % 2. With Covariates
            results = permutation_test_external(diffTrainMatrix, diffTrainBehav, dMatsTr, dBehavTr, ...
                                               diffTestMatrix, dMatsTe, dBehavTe, covars, outputTxt, results);
        else
            fprintf('   (Skipping Block C: Difference files not found)\n');
        end
        
    end
end

%% 3. SAVE RESULTS
fprintf('\nValidation complete. Saving results...\n');
header = {'Training Data', 'Testing  Set', 'r_comb', 'p_comb', 'r_pos', 'p_pos', 'r_neg', 'p_neg', ...
          'mse_comb', 'mse_p_comb', 'mse_pos', 'mse_p_pos', 'mse_neg', 'mse_p_neg'};
      
resultsWithHeader = [header; results];
writecell(resultsWithHeader, outputXlsx);