function [internalTable, crossDrugTable] = build_table_s13(resultRoot)
%BUILD_TABLE_S13 Combine verified 416-, 432-, and 439-node results.
atlases = ["original","no_cerebellum","with_cerebellum"];
atlasLabels = ["Schaefer400+AAL3 (416)", ...
    "Schaefer400+Tian S2 (432)", ...
    "Schaefer400+Tian S2+Buckner7 (439)"];
outcomes = ["GDE","MEQ30","VRS","OBN","BDE","AED"];
internalRows = cell(18,10); outputRow = 0;
for atlasIndex = 1:3
    pathName = fullfile(resultRoot, atlases(atlasIndex), ...
        'internal_validation', 'internal_validation_results.csv');
    values = readtable(pathName, 'TextType','string');
    for outcomeIndex = 1:6
        match = contains(values.Analysis, "LSD_" + outcomes(outcomeIndex));
        row = values(match,:);
        assert(height(row)==1, 'CPM:TableS13:InternalRow', ...
            'Expected one %s row for %s.', outcomes(outcomeIndex), atlases(atlasIndex));
        outputRow = outputRow + 1;
        [cerebellarEdges, cerebellarFraction, withinCerebellum] = ...
            cerebellar_values(resultRoot, atlases(atlasIndex), ...
            outcomes(outcomeIndex));
        internalRows(outputRow,:) = {outcomes(outcomeIndex), atlasLabels(atlasIndex), ...
            row.r_positive, row.p_positive, row.mse_positive, ...
            row.positive_consensus_edges, internal_n(outcomes(outcomeIndex)), ...
            cerebellarEdges, cerebellarFraction, withinCerebellum};
    end
end
internalTable = cell2table(internalRows, 'VariableNames', ...
    {'Outcome','Atlas','R','P','MSE','PositiveConsensusEdges','N', ...
     'CerebellarEdges','CerebellarEdgeFraction','WithinCerebellumEdges'});

crossRows = cell(18,6); outputRow = 0;
for atlasIndex = 1:3
    pathName = fullfile(resultRoot, atlases(atlasIndex), ...
        'cross_drug_loso', 'cross_drug_loso_results.csv');
    values = readtable(pathName, 'TextType','string', ...
        'VariableNamingRule','preserve');
    for outcomeIndex = 1:3
        for testSet = ["other_psych","all_amphs"]
            match = values.("Behavior Scale") == outcomes(outcomeIndex) & ...
                values.("Testing Set") == testSet & values.("P Threshold") == 0.01;
            row = values(match,:);
            assert(height(row)==1, 'CPM:TableS13:CrossDrugRow', ...
                'Expected one %s/%s row for %s.', ...
                outcomes(outcomeIndex), testSet, atlases(atlasIndex));
            outputRow = outputRow + 1;
            crossRows(outputRow,:) = {outcomes(outcomeIndex), atlasLabels(atlasIndex), ...
                testSet, row.r_pos, row.p_pos, row.("N Test")};
        end
    end
end
crossDrugTable = cell2table(crossRows, 'VariableNames', ...
    {'Outcome','Atlas','TestClass','R','P','NTestSessions'});
writetable(internalTable, fullfile(resultRoot, 'table_S13_internal.csv'));
writetable(crossDrugTable, fullfile(resultRoot, 'table_S13_cross_drug.csv'));
end

function value = internal_n(outcome)
if outcome == "MEQ30", value = 48; else, value = 67; end
end

function [edgeCount, fraction, withinCount] = cerebellar_values( ...
        resultRoot, atlasName, outcome)
edgeCount = NaN;
fraction = NaN;
withinCount = NaN;
if atlasName ~= "with_cerebellum"
    return
end
summaryPath = fullfile(resultRoot, atlasName, ...
    'cerebellar_consensus_edge_summary.csv');
summary = readtable(summaryPath, 'TextType','string');
row = summary(summary.Outcome == outcome,:);
assert(height(row)==1, 'CPM:TableS13:CerebellarRow', ...
    'Expected one cerebellar-edge row for %s.', outcome);
edgeCount = row.CerebellarEdges;
fraction = row.CerebellarEdgeFraction;
withinCount = row.WithinCerebellumEdges;
end
