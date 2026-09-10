function [statsA, statsB, audit] = CPM_core_function_cross_drug_loso_reference( ...
        trainEdges, trainBehav, heldoutTrainIdx, ...
        testEdgesA, testBehavA, testEdgesB, testBehavB, ...
        pThreshold, corrType, collectAudit, covarsTrain)
% CPM_CORE_FUNCTION_CROSS_DRUG_LOSO_REFERENCE
%   Conventional CPM inside a subject-wise outer leave-one-out loop.
%   Reference developed during manuscript revision, before caching was added.
%   This is not the consensus-mask external-validation core from v1.0.0.
%
%   Each test subject has one LSD observation in trainEdges and one
%   observation in each of two other-drug test matrices. Before predicting
%   that subject, the matching LSD row is removed from BOTH edge selection
%   and regression fitting. The same fitted fold model predicts both drugs.
%
%   The inputs contain unique upper-triangle edges only. Consequently,
%   network strengths are summed once and are not divided by two.

    if nargin < 10
        collectAudit = false;
    end
    if nargin < 11, covarsTrain = []; end

    validateattributes(trainEdges, {'numeric'}, {'2d', 'real', 'finite'});
    validateattributes(trainBehav, {'numeric'}, {'vector', 'real', 'finite'});
    validateattributes(heldoutTrainIdx, {'numeric'}, ...
        {'vector', 'integer', 'positive'});
    validateattributes(testEdgesA, {'numeric'}, {'2d', 'real', 'finite'});
    validateattributes(testBehavA, {'numeric'}, {'vector', 'real', 'finite'});
    validateattributes(testEdgesB, {'numeric'}, {'2d', 'real', 'finite'});
    validateattributes(testBehavB, {'numeric'}, {'vector', 'real', 'finite'});
    validateattributes(pThreshold, {'numeric'}, ...
        {'scalar', 'real', 'finite', '>', 0, '<', 1});

    corrType = char(validatestring(corrType, {'Pearson', 'Spearman'}));
    trainBehav = trainBehav(:);
    heldoutTrainIdx = heldoutTrainIdx(:);
    testBehavA = testBehavA(:);
    testBehavB = testBehavB(:);

    nEdges = size(trainEdges, 1);
    nTrainTotal = size(trainEdges, 2);
    nTest = numel(heldoutTrainIdx);

    assert(numel(trainBehav) == nTrainTotal, ...
        'CPM:CrossDrugLOSO:TrainAlignment', ...
        'Training edge columns and training behavior are not aligned.');
    assert(size(testEdgesA, 1) == nEdges && size(testEdgesB, 1) == nEdges, ...
        'CPM:CrossDrugLOSO:EdgeCountMismatch', ...
        'Training and testing edge matrices have different edge counts.');
    assert(size(testEdgesA, 2) == nTest && size(testEdgesB, 2) == nTest, ...
        'CPM:CrossDrugLOSO:TestAlignment', ...
        'Each test drug must have one column per held-out subject.');
    assert(numel(testBehavA) == nTest && numel(testBehavB) == nTest, ...
        'CPM:CrossDrugLOSO:TestBehaviorAlignment', ...
        'Testing edges and behaviors are not aligned.');
    assert(numel(unique(heldoutTrainIdx)) == nTest, ...
        'CPM:CrossDrugLOSO:DuplicateHeldoutSubject', ...
        'Each target-cohort subject must map to one unique LSD row.');
    assert(all(heldoutTrainIdx <= nTrainTotal), ...
        'CPM:CrossDrugLOSO:HeldoutIndexOutOfRange', ...
        'A held-out LSD index lies outside the eligible training set.');
    assert(std(testBehavA) > 0 && std(testBehavB) > 0, ...
        'CPM:CrossDrugLOSO:ConstantTestBehavior', ...
        'A test-drug behavior vector is constant.');

    predAComb = zeros(nTest, 1);
    predAPos = zeros(nTest, 1);
    predANeg = zeros(nTest, 1);
    predBComb = zeros(nTest, 1);
    predBPos = zeros(nTest, 1);
    predBNeg = zeros(nTest, 1);
    baselinePred = zeros(nTest, 1);

    if collectAudit
        posSelectionCount = zeros(nEdges, 1, 'uint16');
        negSelectionCount = zeros(nEdges, 1, 'uint16');
        posEdgeCountByFold = zeros(nTest, 1);
        negEdgeCountByFold = zeros(nTest, 1);
    end

    for fold = 1:nTest
        keepTrain = true(nTrainTotal, 1);
        keepTrain(heldoutTrainIdx(fold)) = false;

        foldEdges = trainEdges(:, keepTrain);
        foldBehav = trainBehav(keepTrain);

        % Direct feature selection on all n-1 eligible LSD subjects.
        if isempty(covarsTrain)
            [edgeR, edgeP] = corr(foldEdges', foldBehav, 'Type', corrType);
        else
            [edgeR, edgeP] = partialcorr(foldEdges', foldBehav, ...
                covarsTrain(keepTrain,:), 'Type', corrType);
        end
        posMask = edgeR > 0 & edgeP < pThreshold;
        negMask = edgeR < 0 & edgeP < pThreshold;

        trainSumPos = sum(foldEdges(posMask, :), 1)';
        trainSumNeg = sum(foldEdges(negMask, :), 1)';

        % pinv gives the ordinary least-squares solution and remains defined
        % if a fold has an empty network or collinear network strengths.
        bComb = pinv([trainSumPos, trainSumNeg, ones(sum(keepTrain), 1)]) ...
            * foldBehav;
        bPos = pinv([trainSumPos, ones(sum(keepTrain), 1)]) * foldBehav;
        bNeg = pinv([trainSumNeg, ones(sum(keepTrain), 1)]) * foldBehav;

        testASumPos = sum(testEdgesA(posMask, fold), 1);
        testASumNeg = sum(testEdgesA(negMask, fold), 1);
        testBSumPos = sum(testEdgesB(posMask, fold), 1);
        testBSumNeg = sum(testEdgesB(negMask, fold), 1);

        predAComb(fold) = bComb(1) * testASumPos + ...
            bComb(2) * testASumNeg + bComb(3);
        predAPos(fold) = bPos(1) * testASumPos + bPos(2);
        predANeg(fold) = bNeg(1) * testASumNeg + bNeg(2);

        predBComb(fold) = bComb(1) * testBSumPos + ...
            bComb(2) * testBSumNeg + bComb(3);
        predBPos(fold) = bPos(1) * testBSumPos + bPos(2);
        predBNeg(fold) = bNeg(1) * testBSumNeg + bNeg(2);

        baselinePred(fold) = mean(foldBehav);

        if collectAudit
            posSelectionCount = posSelectionCount + uint16(posMask);
            negSelectionCount = negSelectionCount + uint16(negMask);
            posEdgeCountByFold(fold) = nnz(posMask);
            negEdgeCountByFold(fold) = nnz(negMask);
        end
    end

    statsA = evaluate_predictions( ...
        predAComb, predAPos, predANeg, baselinePred, testBehavA);
    statsB = evaluate_predictions( ...
        predBComb, predBPos, predBNeg, baselinePred, testBehavB);

    assert(all(isfinite([statsA.r, statsA.mse, statsA.baseline_mse, ...
                         statsB.r, statsB.mse, statsB.baseline_mse])), ...
        'CPM:CrossDrugLOSO:NonfinitePerformance', ...
        'The LOSO procedure produced a nonfinite performance statistic.');

    if collectAudit
        audit.pos_selection_count = posSelectionCount;
        audit.neg_selection_count = negSelectionCount;
        audit.pos_edge_count_by_fold = posEdgeCountByFold;
        audit.neg_edge_count_by_fold = negEdgeCountByFold;
        audit.n_train_per_fold = nTrainTotal - 1;
        audit.n_test = nTest;
    else
        audit = struct();
    end
end

function stats = evaluate_predictions(predComb, predPos, predNeg, ...
        baselinePred, observed)
    stats.r = [corr(predComb, observed, 'Type', 'Pearson'), ...
               corr(predPos, observed, 'Type', 'Pearson'), ...
               corr(predNeg, observed, 'Type', 'Pearson')];
    stats.mse = [mean((observed - predComb).^2), ...
                 mean((observed - predPos).^2), ...
                 mean((observed - predNeg).^2)];
    stats.baseline_mse = mean((observed - baselinePred).^2);
    stats.pred_comb = predComb;
    stats.pred_pos = predPos;
    stats.pred_neg = predNeg;
    stats.baseline_pred = baselinePred;
    stats.observed = observed;
end
