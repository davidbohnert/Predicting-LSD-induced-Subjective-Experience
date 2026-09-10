function test_loso_equivalence
%TEST_LOSO_EQUIVALENCE Compare readable and cached subject-wise LOSO cores.
testDir = fileparts(mfilename('fullpath'));
repositoryRoot = fileparts(testDir);
addpath(fullfile(repositoryRoot, 'tests', 'reference'), ...
    fullfile(repositoryRoot, 'scripts', 'common'), ...
    fullfile(repositoryRoot, 'pipelines', 'cross_drug_loso'));

rng(902, 'twister');
nNodes = 9; nTrain = 15; nTest = 5;
trainMats = symmetric_mats(nNodes, nTrain);
testMatsA = symmetric_mats(nNodes, nTest);
testMatsB = symmetric_mats(nNodes, nTest);
trainBehavior = randn(nTrain, 1);
testBehaviorA = randn(nTest, 1);
testBehaviorB = randn(nTest, 1);
heldout = (6:10)';
cohorts = [repmat("5HT",5,1); repmat("LAM",5,1); repmat("LPM",5,1)];
globalRows = (1:nTrain)';
threshold = 0.20;
covars = [double((1:nTrain)' > 7), randn(nTrain,1)];

upper = find(triu(true(nNodes), 1));
trainEdges = reshape(trainMats, nNodes*nNodes, nTrain);
trainEdges = trainEdges(upper, :);
testEdgesA = reshape(testMatsA, nNodes*nNodes, nTest);
testEdgesA = testEdgesA(upper, :);
testEdgesB = reshape(testMatsB, nNodes*nNodes, nTest);
testEdgesB = testEdgesB(upper, :);
[referenceA, referenceB] = CPM_core_function_cross_drug_loso_reference( ...
    trainEdges, trainBehavior, heldout, testEdgesA, testBehaviorA, ...
    testEdgesB, testBehaviorB, threshold, 'Pearson', false);

optimized = permutation_test_cross_drug_loso( ...
    trainMats, trainBehavior, [], globalRows, cohorts, heldout, "LAM", ...
    testMatsA, testBehaviorA, testMatsB, testBehaviorB, threshold, ...
    'Pearson', 2, false);
assert(max(abs(optimized.drug_a.r - referenceA.r)) < 1e-10);
assert(max(abs(optimized.drug_b.r - referenceB.r)) < 1e-10);
assert(max(abs(optimized.drug_a.mse - referenceA.mse)) < 1e-10);
assert(max(abs(optimized.drug_b.mse - referenceB.mse)) < 1e-10);
fprintf('LOSO optimized/reference equivalence passed.\n');

[referenceAAdjusted, referenceBAdjusted] = ...
    CPM_core_function_cross_drug_loso_reference( ...
        trainEdges, trainBehavior, heldout, testEdgesA, testBehaviorA, ...
        testEdgesB, testBehaviorB, threshold, 'Pearson', false, covars);
optimizedAdjusted = permutation_test_cross_drug_loso( ...
    trainMats, trainBehavior, covars, globalRows, cohorts, heldout, "LAM", ...
    testMatsA, testBehaviorA, testMatsB, testBehaviorB, threshold, ...
    'Pearson', 2, false);
assert(max(abs(optimizedAdjusted.drug_a.r-referenceAAdjusted.r)) < 1e-10);
assert(max(abs(optimizedAdjusted.drug_b.r-referenceBAdjusted.r)) < 1e-10);
fprintf('Adjusted LOSO optimized/reference equivalence passed.\n');
end

function mats = symmetric_mats(nNodes, nSubjects)
mats = zeros(nNodes, nNodes, nSubjects);
for subject = 1:nSubjects
    raw = randn(nNodes);
    current = (raw + raw') / 2;
    current(1:nNodes+1:end) = 0;
    mats(:,:,subject) = current;
end
end
